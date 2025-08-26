from chimerax.ui.gui import MainToolWindow
from chimerax.core.tools import ToolInstance
from chimerax.atomic import AtomicStructure, selected_atoms
from chimerax.core.commands import run
from Qt.QtWidgets import QWidget, QVBoxLayout, QGridLayout, QPushButton, QLabel, QComboBox, QCheckBox, QFileDialog, QLineEdit
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import mdtraj as md
import numpy as np

class AmberDistTool(ToolInstance):
    SESSION_ENDURING = True
    SESSION_SAVE = True

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)

        # Main container
        self.tool_window = MainToolWindow(self)

        layout = QVBoxLayout()
        entry_layout = QGridLayout()

        # File selector
        self.traj_label = QLabel("Trajectory File")
        self.traj_name = QLineEdit()
        self.traj_name.setMinimumWidth(50)
        self.traj_name.setPlaceholderText("select file")
        btn_traj_file = QPushButton("Choose trajectory (.nc, .dcd, etc.)")
        btn_traj_file.clicked.connect(self.choose_traj_file)
        entry_layout.addWidget(self.traj_label, 0, 0)
        entry_layout.addWidget(self.traj_name, 0, 1)
        entry_layout.addWidget(btn_traj_file, 0, 2)

        self.parm_label = QLabel("Topology File")
        self.parm_name = QLineEdit()
        self.parm_name.setMinimumWidth(50)
        self.parm_name.setPlaceholderText("select file")
        btn_parm_file = QPushButton("Choose topology (.parm7, .pdb, etc.")
        btn_parm_file.clicked.connect(self.choose_parm_file)
        entry_layout.addWidget(self.parm_label, 1, 0)
        entry_layout.addWidget(self.parm_name, 1, 1)
        entry_layout.addWidget(btn_parm_file, 1, 2)
       
        self.atom1_label = QLabel("Atom 1")
        self.atom1_entry = QLineEdit()
        self.atom1_entry.setPlaceholderText("required")
        self.atom1_label.hide()
        self.atom1_entry.hide()
        entry_layout.addWidget(self.atom1_label, 2, 0)
        entry_layout.addWidget(self.atom1_entry, 3, 0)
        self.atom2_label = QLabel("Atom 2")
        self.atom2_entry = QLineEdit()
        self.atom2_entry.setPlaceholderText("required")
        self.atom2_label.hide()
        self.atom2_entry.hide()
        entry_layout.addWidget(self.atom2_label, 2, 1)
        entry_layout.addWidget(self.atom2_entry, 3, 1)
        self.atom3_label = QLabel("Atom 3")
        self.atom3_entry = QLineEdit()
        self.atom3_entry.setPlaceholderText("for angles")
        self.atom3_label.hide()
        self.atom3_entry.hide()
        entry_layout.addWidget(self.atom3_label, 2, 2)
        entry_layout.addWidget(self.atom3_entry, 3, 2)
        
        self.sel_label = QLabel("Selection")
        self.sel_entry = QLineEdit()
        self.sel_entry.setPlaceholderText("default is protein")
        self.sel_label.hide()
        self.sel_entry.hide()
        entry_layout.addWidget(self.sel_label, 2, 0)
        entry_layout.addWidget(self.sel_entry, 2, 1)

        layout.addLayout(entry_layout)

        # Analysis type
        self.combo = QComboBox()
        self.combo.addItems(["--Select--", "Geometries", "Hydrogen Bonding"])
        layout.addWidget(QLabel("Analysis:"))
        layout.addWidget(self.combo)
        self.combo.currentTextChanged.connect(self.toggle_combo)

        self.hbond_combo = QComboBox()
        self.hbond_combo.addItems(["Across Trajectory", "Contact Map"])
        self.hbond_combo_lab = QLabel("H-Bond Analysis Type:")
        self.hbond_combo.hide()
        self.hbond_combo_lab.hide()
        layout.addWidget(self.hbond_combo_lab)
        layout.addWidget(self.hbond_combo)
        
        
        # Compute button
        btn_compute = QPushButton("Compute")
        btn_compute.clicked.connect(self.compute)
        layout.addWidget(btn_compute)

        # Matplotlib canvas
        self.figure = Figure(figsize=(5,4), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)

        # Attach layout to tool window
        self.tool_window.ui_area.setLayout(layout)

        self.tool_window.manage(placement="side")
        # Trajectory file placeholder
        self.trajfile = None
        
    def toggle_combo(self, stat):
        #stat = self.combo.currentText()
        if stat == "Geometries":
            self.atom1_label.show()
            self.atom2_label.show()
            self.atom3_label.show()
            self.atom1_entry.show()
            self.atom2_entry.show()
            self.atom3_entry.show()
            self.hbond_combo_lab.hide()
            self.hbond_combo.hide()
        elif stat == "Hydrogen Bonding":
            self.atom1_label.hide()
            self.atom2_label.hide()
            self.atom3_label.hide()
            self.atom1_entry.hide()
            self.atom2_entry.hide()
            self.atom3_entry.hide()
            self.sel_label.show()
            self.sel_entry.show()
            self.hbond_combo.show()
            self.hbond_combo_lab.show()

    def choose_traj_file(self):
        fname, _ = QFileDialog.getOpenFileName(None, "Open Trajectory", "", "Trajectory Files (*.nc *.dcd *.xtc);;All Files (*)")
        if fname:
            self.trajfile = fname
            self.traj_name.setText(fname)

    def choose_parm_file(self):
        fname, _ = QFileDialog.getOpenFileName(None, "Open Topology", "", "Topology Files (*.parm7 *.pdb);;All Files (*)")
        if fname:
            self.top = fname
            if fname.split(".")[-1] == "pdb":
                run(self.session, f"open {fname}")
            else:
                import tempfile
                import parmed as pmd
                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as temp_file:
                    temp_filepath = temp_file.name
                amber = pmd.load_file(fname, self.trajfile)
                amber.coordinates = amber.coordinates[:,:]
                amber.save(temp_filepath, format="pdb", charmm=True, overwrite=True, renumber=True)
                run(self.session, f"open {temp_filepath}")
            self.parm_name.setText(fname)
   
    def compute(self):
        if not self.trajfile:
            self.session.logger.warning("No trajectory file selected.")
            return
        

        type_anal = self.combo.currentText()

        traj = md.load(self.trajfile, top=self.top)

        if type_anal == "Geometries":
            self.sel_atoms = selected_atoms(self.session)
            if len(self.sel_atoms) > 0:
                self.indices = [a.coord_index for a in self.sel_atoms]
                self.atom1_entry.setText(str(self.indices[0]))
                self.atom2_entry.setText(str(self.indices[1]))
                if len(self.sel_atoms) > 2:
                    self.atom3_entry.setText(str(self.indices[2]))
                if len(self.sel_atoms) > 3:
                    self.session.logger.warning("Only select a maximum of 3 atoms")
            
            if len(self.sel_atoms) == 2:
                self.figure.clear()
                ax = self.figure.add_subplot(111)
                values = md.compute_distances(traj, [self.indices])
                values_a = values * 10
                ax.plot(values_a)
                ax.set_xlabel("Frame")
                ax.set_ylabel("Distance (Å)")
                ax.set_title(f"Distances between {self.indices[0]} and {self.indices[1]} over time")

            elif len(self.sel_atoms) == 3:
                self.figure.clear()
                ax = self.figure.add_subplot(111)
                values = md.compute_angles(traj, [self.indices])
                angles_deg = np.degrees(values)
                ax.plot(angles_deg)
                ax.set_xlabel("Frame")
                ax.set_ylabel("Angle (°)")
                ax.set_title(f"Angels between {self.indices} over time")
        elif type_anal == "Hydrogen Bonding":
            hbond_stat = self.hbond_combo.currentText()
            sel_atoms = selected_atoms(self.session)

            if len(sel_atoms) > 0:
                sel = [a.coord_index for a in sel_atoms]
                self.sel_entry.setText("Using selected")
            else:
                selection = self.sel_entry.text().strip()
                if selection and selection not in ("protein", "backbone"):
                    atoms = f"name {selection}"
                else:
                    atoms = selection or "protein"
                sel = traj.topology.select(atoms)

            if hbond_stat == "Across Trajectory":
                self.figure.clear()
                ax = self.figure.add_subplot(111)
                hbonds = md.baker_hubbard(traj, freq=0.0)  # geometric criterion, all hbonds
                hb_per_frame_sel = []
                for i in range(traj.n_frames):
                    hbonds = md.baker_hubbard(traj[i], freq=0.3)
                    # keep only hbonds where both donor/acceptor are in selection
                    hbonds_sel = [hb for hb in hbonds if hb[0] in sel or hb[2] in sel]
                    hb_per_frame_sel.append(len(hbonds_sel))
                hb_per_frame_sel = np.array(hb_per_frame_sel)
                ax.plot(hb_per_frame_sel)
                ax.set_xlabel("Frame")
                ax.set_ylabel("# fo H-Bonds")
                ax.set_title(f"# of Hydrogen Bonds over time")
            elif hbond_stat == "Contact Map":
                self.figure.clear()
                ax = self.figure.add_subplot(111)
                hbonds = md.baker_hubbard(traj, freq=0.0)  # list of (donor, hydrogen, acceptor)

                def hbond_occupancy(traj, hbonds, distance_cutoff=0.35, angle_cutoff=2.0944):
                    donor_hydrogen = hbonds[:, [0, 1]]
                    hydrogen_acceptor = hbonds[:, [1, 2]]

                    dists = md.compute_distances(traj, hydrogen_acceptor)
                    angles = np.degrees(md.compute_angles(traj, hbonds))

                    present = (dists < distance_cutoff) & (angles > angle_cutoff)
                    return present
                hb_array = hbond_occupancy(traj, hbonds)  # shape (n_frames, n_hbonds)
                occupancy = hb_array.mean(axis=0)  # per bond
                
                sel_residues = sorted({traj.topology.atom(i).residue.index for i in sel})
                res_map = {res: i for i, res in enumerate(sel_residues)}
                sel_hbonds = [
                       (d, h, a) for (d, h, a), frac in zip(hbonds, occupancy)
                        if d in sel and a in sel
                        ]

                sel_occupancy = [
                        frac for (d, h, a), frac in zip(hbonds, occupancy)
                        if d in sel and a in sel
                        ]
                
                atom_to_res = np.array([atom.residue.index for atom in traj.topology.atoms])
                donors = [hb[0] for hb in hbonds]
                acceptors = [hb[2] for hb in hbonds]
                n_sel_res = len(sel_residues)
                sel_contact_map = np.zeros((n_sel_res, n_sel_res))
                #n_res = len(sel) #traj.topology.n_residues
                #res_contact_map = np.zeros((n_res, n_res))
                #for frac, d, a in zip(occupancy, donors, acceptors):
                #for frac, (d, h, a) in zip(sel_occupancy, sel_hbonds):
                #    res_d = atom_to_res[d]
                #    res_a = atom_to_res[a]
                #    res_contact_map[res_d, res_a] = max(res_contact_map[res_d, res_a], frac)
                #    res_contact_map[res_a, res_d] = res_contact_map[res_d, res_a] 
                for frac, (d, h, a) in zip(sel_occupancy, sel_hbonds):
                    res_d = atom_to_res[d]
                    res_a = atom_to_res[a]
                    if res_d in res_map and res_a in res_map:
                        i = res_map[res_d]
                        j = res_map[res_a]
                        sel_contact_map[i, j] = max(sel_contact_map[i, j], frac)
                        sel_contact_map[j, i] = sel_contact_map[i, j]
                im = ax.imshow(sel_contact_map, cmap="viridis", origin="lower", vmin=0, vmax=1)
                self.figure.colorbar(im, ax=ax, label="Occupancy fraction")
                ax.set_xlabel("Acceptor Residue index")
                ax.set_ylabel("Donor Residue index")
                ax.set_title("Hydrogen Bond Contact Map")

        self.canvas.draw()

