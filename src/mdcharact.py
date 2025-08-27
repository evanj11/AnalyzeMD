from chimerax.ui.gui import MainToolWindow
from chimerax.core.tools import ToolInstance
from chimerax.atomic import AtomicStructure, selected_atoms
from chimerax.core.commands import run
from Qt.QtWidgets import QWidget, QVBoxLayout, QGridLayout, QPushButton, QLabel, QComboBox, QCheckBox, QFileDialog, QLineEdit
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import mdtraj as md
import numpy as np

class AmberMDTool(ToolInstance):
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
        layout.addLayout(entry_layout)

        param_layout = QGridLayout()

        self.sel_label = QLabel("Selection:")
        self.sel_entry = QLineEdit()
        self.sel_entry.setPlaceholderText("default protein")
        param_layout.addWidget(self.sel_label, 0, 0)
        param_layout.addWidget(self.sel_entry, 0, 1)

        #self.play_btn = QPushButton("Play")
        #self.play_btn.clicked.connect(self.play)
        #entry_layout.addWidget(self.play_btn, 2, 2)

        self.stride = QLineEdit()
        self.stride.setPlaceholderText("default is 1")
        self.stride_label = QLabel("Stride:")
        param_layout.addWidget(self.stride_label, 0, 2)
        param_layout.addWidget(self.stride, 0, 3)
        layout.addLayout(param_layout)

        # Analysis type
        anal_layout = QGridLayout()
        self.combo_label = QLabel("Analysis:")
        self.combo = QComboBox()
        self.combo.addItems(["RMSD", "RMSF", "RadG"])
        anal_layout.addWidget(self.combo_label, 3, 0)
        anal_layout.addWidget(self.combo, 3, 1)
        layout.addLayout(anal_layout)

        #self.pause_btn = QPushButton("Pause")
        #self.pause_btn.clicked.connect(self.pause)
        #entry_layout.addWidget(self.pause_btn, 3, 2)


        # Free energy option
        self.free_energy = QCheckBox("Calculate Free Energy")
        self.free_energy.clicked.connect(self.free_ene_toggle)
        layout.addWidget(self.free_energy)
        
        temp_layout = QGridLayout()
        self.temp_label = QLabel("Temperature (if known) for FE")
        self.temp_entry = QLineEdit()
        self.temp_label.hide()
        self.temp_entry.hide()
        temp_layout.addWidget(self.temp_label, 0, 0)
        temp_layout.addWidget(self.temp_entry, 0, 1)
        layout.addLayout(temp_layout)

        self.pairwise = QCheckBox("Pairwise Graph?")
        self.pairwise.clicked.connect(self.open_pairwise)
        layout.addWidget(self.pairwise)
        
        self.axes_layout = QGridLayout()
        self.axis1_label = QLabel("Axes for Pairwise Graph:")
        self.axis1 = QComboBox()
        self.axis1.addItems(["RMSD", "PCA", "DCCM"])
        self.axis1.hide()
        self.axis1_label.hide()
        self.axes_layout.addWidget(self.axis1_label, 0, 0)
        self.axes_layout.addWidget(self.axis1, 0, 1)
        layout.addLayout(self.axes_layout)
        
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


    def free_ene_toggle(self):
        if self.free_energy.isChecked():
            self.temp_entry.show()
            self.temp_label.show()
            self.combo_label.hide()
            self.combo.hide()
            self.axis1.hide()
            self.axis1_label.hide()
            self.pairwise.setChecked(False)
        else:
            self.temp_entry.hide()
            self.temp_label.hide()
            self.combo.show()
            self.combo_label.show()

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

    def play(self):
        run(self.session, f"open {self.trajfile}")

    def open_pairwise(self):
        if self.pairwise.isChecked():
            self.axis1.show()
            self.axis1_label.show()
            self.combo.hide()
            self.combo_label.hide()
            self.temp_entry.hide()
            self.temp_label.hide()
            self.free_energy.setChecked(False)
        else:
            self.axis1.hide()
            self.axis1_label.hide()
            self.combo.show()
            self.combo_label.show()

    def compute(self):
        if not self.trajfile:
            self.session.logger.warning("No trajectory file selected.")
            return

        stat = self.combo.currentText()
        fe_flag = self.free_energy.isChecked()
        pairwise_flag = self.pairwise.isChecked()

        if self.stride.text().strip():
            stride = int(self.stride.text().strip())
        else:
            stride = 1

        traj = md.load(self.trajfile, stride=stride, top=self.top)
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
     
        if stat == "RMSD":
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            values = md.rmsd(traj, traj, frame=0, atom_indices=indices)
            values_a = values * 10
            ax.plot(values_a)
            ax.set_xlabel("Frame")
            ax.set_ylabel("RMSD (Å)")
            ax.set_title("RMSD over time")

        elif stat == "RMSF":
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            values = md.rmsf(traj, traj, 0, atom_indices=indices)
            values_a = values * 10
            ax.bar(np.arange(len(values_a)), values_a)
            ax.set_xlabel("Atom index")
            ax.set_ylabel("RMSF (Å)")
            ax.set_title("RMSF per atom")

        elif stat == "RadG":
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            traj_slice = traj.atom_slice(indices)
            values = md.compute_rg(traj_slice)
            values_a = values * 10
            ax.plot(values_a)
            ax.set_xlabel("Frame")
            ax.set_ylabel("Radius of gyration (Å)")
            ax.set_title("Radius of gyration")
        
        self.canvas.draw()

        # Free energy surface
        if fe_flag:
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            
            traj_slice = traj.atom_slice(indices)
            rmsd = md.rmsd(traj_slice, traj_slice, 0)
            rmsd = rmsd * 10
            rg = md.compute_rg(traj_slice)
            rg = rg * 10

            H, xedges, yedges = np.histogram2d(rmsd, rg, bins=50)
            if self.temp_entry.text().strip():
                T = float(self.temp_entry.text().strip())
                F = -0.0083145 * T * np.log(H + 1e-12)
            
                F[F > 25] = np.nan
                im = ax.contourf(xedges[:-1], yedges[:-1], F.T, levels=50, cmap="seismic")
                ax.set_xlabel("RMSD (Å)")
                ax.set_ylabel("Rg (Å)")
                ax.set_title("Free Energy Landscape")
                self.figure.colorbar(im, ax=ax, label="ΔG (kcal/mol)")
                self.canvas.draw()
            else:
                F = -np.log(H + 1e-12)
            
                F[F > 25] = np.nan
                im = ax.contourf(xedges[:-1], yedges[:-1], F.T, levels=50, cmap="seismic")
                ax.set_xlabel("RMSD (Å)")
                ax.set_ylabel("Rg (Å)")
                ax.set_title("Free Energy Landscape")
                self.figure.colorbar(im, ax=ax, label="ΔG (kBT)")
                self.canvas.draw()

        if pairwise_flag:
            self.figure.clear()
            ax = self.figure.add_subplot(111)

            if self.axis1.currentText() == "RMSD":
                traj_slice = traj.atom_slice(indices)
                n_frames = traj.n_frames
                rmsd_matrix = np.zeros((n_frames, n_frames))

                for i in range(n_frames):
                    rmsd_matrix[i, :] = md.rmsd(traj_slice, traj_slice, i) * 10

                im = ax.imshow(rmsd_matrix, cmap="viridis", origin="lower")
                self.figure.colorbar(im, ax=ax, label="RMSD (Å)")
                ax.set_xlabel("Frame")
                ax.set_ylabel("Frame")
                ax.set_title("Pairwise RMSD Matrix")
                self.canvas.draw()

            elif self.axis1.currentText() == "PCA":
                from sklearn.decomposition import PCA
                atom_indices = indices

                X = traj.xyz[:, atom_indices, :].reshape(traj.n_frames, -1)
                X -= X.mean(axis=0)

                pca = PCA(n_components=2)
                X_pca = pca.fit_transform(X)
                im = ax.scatter(X_pca[:,0], X_pca[:,1], s=10, alpha=0.7, c=range(traj.n_frames), cmap="viridis")
                ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
                ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
                self.figure.colorbar(im, ax=ax, label="Frame index")
                ax.set_title("PCA of MD trajectory")
                def on_pick(event):
                    ind = event.ind[0]   # index of picked frame
                    frame_idx = ind
                    # Load that frame into ChimeraX viewer
                    traj[frame_idx].save_pdb("picked_frame.pdb")
                    run(self.session, f"open picked_frame.pdb")

                im.set_picker(5)   # tolerance in points
                self.canvas.mpl_connect('pick_event', on_pick)
                self.canvas.draw()
            
            elif self.axis1.currentText() == "DCCM":
                atom_indices = indices

                X = traj.xyz[:, atom_indices, :]
                X -= X.mean(axis=0)
                n_atoms = X.shape[1]
                dccm = np.zeros((n_atoms, n_atoms))
                for i in range(n_atoms):
                    for j in range(n_atoms):
                        num = np.sum(np.mean(X[:,i,:] * X[:,j,:], axis=0))   # dot product
                        den = np.sqrt(
                            np.sum(np.mean(X[:,i,:] * X[:,i,:], axis=0)) *
                            np.sum(np.mean(X[:,j,:] * X[:,j,:], axis=0))
                            )
                        dccm[i,j] = num / den
                im = ax.imshow(dccm, cmap="RdBu_r", vmin=-1, vmax=1, origin="lower")
                self.figure.colorbar(im, ax=ax, label="Cross-Correlation")
                ax.set_xlabel("Residue Index")
                ax.set_ylabel("Residue Index")
                ax.set_title("DCCM Matrix")
                self.canvas.draw()

            else:
                self.session.logger.warning("Please select both axes the same")
