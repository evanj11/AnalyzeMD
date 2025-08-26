from chimerax.ui.gui import MainToolWindow
from chimerax.core.tools import ToolInstance
from chimerax.atomic import AtomicStructure, selected_atoms
from chimerax.core.commands import run
from Qt.QtWidgets import QWidget, QVBoxLayout, QGridLayout, QPushButton, QLabel, QComboBox, QCheckBox, QFileDialog, QLineEdit
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import mdtraj as md
import numpy as np
from deeptime.decomposition import TICA, VAMP
from deeptime.clustering import KMeans
from deeptime.markov.msm import MaximumLikelihoodMSM
from deeptime.markov import pcca

class AmberMSMTool(ToolInstance):
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

        self.sel_label = QLabel("Selection")
        self.sel_entry = QLineEdit()
        self.sel_entry.setPlaceholderText("default protein, choose CA, backbone, etc.")
        entry_layout.addWidget(self.sel_label, 2, 0)
        entry_layout.addWidget(self.sel_entry, 2, 1)
        
        # Analysis type
        self.combo_label = QLabel("Features:")
        self.combo = QComboBox()
        self.combo.addItems(["RMSD-RadG", "Phi-Psi"])
        entry_layout.addWidget(self.combo_label, 3, 0)
        entry_layout.addWidget(self.combo, 3, 1)

        self.number_micstates_label = QLabel("Number of Microstates to Compute")
        self.number_micstates = QLineEdit()
        self.number_micstates.setPlaceholderText("default is 100, must be > #frames")
        entry_layout.addWidget(self.number_micstates_label, 4, 0)
        entry_layout.addWidget(self.number_micstates, 4, 1)

        # Free energy option
        self.number_macstates_label = QLabel("Number of Macrostates to Compute")
        self.number_macstates = QLineEdit()
        entry_layout.addWidget(self.number_macstates_label, 5, 0)
        entry_layout.addWidget(self.number_macstates, 5, 1)
        
        layout.addLayout(entry_layout)
        
        btn_compute = QPushButton("Compute")
        btn_compute.clicked.connect(self.compute)
        layout.addWidget(btn_compute)

        # Matplotlib canvas
        self.figure = Figure(figsize=(5,4), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        
        self.show_models = QPushButton("Show Macrostate Models")
        self.show_models.clicked.connect(self.show_mods)
        layout.addWidget(self.show_models)

        # Attach layout to tool window
        self.tool_window.ui_area.setLayout(layout)

        self.tool_window.manage(placement="side")
        # Trajectory file placeholder
        self.trajfile = None

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

        stat = self.combo.currentText()

        traj = md.load(self.trajfile, top=self.top)
        sel_atoms = selected_atoms(self.session)

        if len(sel_atoms) > 0:
            indices = [a.coord_index for a in sel_atoms]
        else:
            selection = self.sel_entry.text().strip()
            if selection and selection not in ("protein", "backbone"):
                atoms = f"name {selection}"
            else:
                atoms = selection or "protein"
            indices = traj.topology.select(atoms)

        traj.center_coordinates()
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        traj_slice = traj.atom_slice(indices)

        if stat == "Phi-Psi":
            phi_indices, phi_angles = md.compute_phi(traj_slice)
            psi_indices, psi_angles = md.compute_psi(traj_slice)
            features = np.hstack([phi_angles, psi_angles])
        elif stat == "RMSD-RadG":
            rmsd = md.rmsd(traj_slice, traj_slice, 0)
            radg = md.compute_rg(traj_slice)
            features = np.column_stack([rmsd, radg])  
        
        self.session.logger.warning(f"{features}")

        tica = TICA(lagtime=5, dim=15)   # adjust lagtime
        tica_data = tica.fit_transform(features)

#        vamp = VAMP(lagtime=5)
#        tica_data = vamp.fit_transform(features)

        if self.number_micstates.text().strip():
            n_clusters = int(self.number_micstates.text().strip())
        else:
            n_clusers = traj.n_frames / 1000

        kmeans = KMeans(n_clusters=n_clusters, max_iter=500)
        clustering = KMeans(n_clusters=n_clusters).fit(tica_data).fetch_model()
        dtrajs = clustering.transform(tica_data)
        
        msm = MaximumLikelihoodMSM().fit(dtrajs, lagtime=1).fetch_model()

        m = int(self.number_macstates.text().strip())
        
        T = msm.transition_matrix              # (n_microstates x n_microstates)
        pi = msm.stationary_distribution       # stationary distribution

        pcca = msm.pcca(n_metastable_sets=m)
        
        memberships = pcca.memberships         # shape (n_microstates, m)

        state_to_macro = memberships.argmax(axis=1)

        macro_of_frame = state_to_macro[dtrajs]
        
        mu = pi
        macro_sets = [np.where(state_to_macro == k)[0] for k in range(m)]

        rep_micro_per_macro = [S[np.argmax(mu[S])] for S in macro_sets]

        rep_frame_per_macro = []
        for k, s_star in enumerate(rep_micro_per_macro):
            frames_in_s = np.where(dtrajs == s_star)[0]
            if len(frames_in_s) == 0:
                rep_frame_per_macro.append(None)
            else:
                rep_frame_per_macro.append(frames_in_s[len(frames_in_s)//2])

        for k, fidx in enumerate(rep_frame_per_macro):
            if fidx is not None:
                traj[fidx].save_pdb(f"macro_{k}_rep.pdb")

        macro_trajs = [state_to_macro[dtraj] for dtraj in dtrajs]
        projected_data = np.vstack(tica_data)  # Your original 2D data
        macro_of_frame = np.array(state_to_macro[dtrajs])  # (n_frames,)
        macro_data = [projected_data[macro_of_frame == i] for i in range(m)]
        
        A = [0]   # index or set of macrostates
        B = [m]

        #tpt_result = tpt(msm.transition_matrix, msm.stationary_distribution, A, B)
        flux = msm.reactive_flux(A, B)
        paths, capacities = flux.pathways(fraction=.8, maxiter=1000)
        best_path = paths[0]       # list of state indices (microstates or macrostates depending on definition)
        best_capacity = capacities[0]
        state_coords = []
        for s in best_path:
            frames_in_s = np.where(dtrajs == s)[0]
            if len(frames_in_s) > 0:
                mean_coord = tica_data[frames_in_s].mean(axis=0)[:2]  # TC1, TC2
                state_coords.append(mean_coord)

        colors = ['red', 'blue', 'green', 'orange', 'magenta', 'purple', 'deepskyblue', 'darkgray', 'yellow']
        for i, data in enumerate(macro_data):
            im = ax.scatter(data[:, 0], data[:, 1], s=5, color=colors[i], label=f'State {i}', alpha=0.5)    
        
        state_coords = np.array(state_coords)
        ax.plot(state_coords[:,0], state_coords[:,1], '-x', color='black', lw=1, markersize=4, label=f"Capacity {best_capacity:.2e}")

        ax.set_xlabel("TIC 1")
        ax.set_ylabel("TIC 2")
        ax.set_title("Macrostates (PCCA+) in TICA space")
        ax.legend()
        self.canvas.draw()

    def show_mods(self):
        colors = ['indianred', 'royalblue', 'forestgreen', 'orange', 'orchid', 'darkorchid', 'lightskyblue', 'gray', 'khaki']
        m = int(self.number_macstates.text().strip())
        run(self.session, f"hide #1 model")
        mods = []
        for i in range(m):
            run(self.session, f"open macro_{i}_rep.pdb")
            run(self.session, f"color #{i+2} {colors[i]}")
            mods.append(i+2)
        models = ",".join(map(str, mods))
        run(self.session, f"matchmaker #{models} to #1")
        run(self.session, f"morph #{models} frames 40")
        run(self.session, f"hide #{i+3} model")
        run(self.session, f"show #{models} model")
