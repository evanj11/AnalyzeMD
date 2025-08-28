from chimerax.ui.gui import MainToolWindow
from chimerax.core.tools import ToolInstance
from chimerax.atomic import AtomicStructure, selected_atoms
from chimerax.core.commands import run
from chimerax.markers import MarkerSet, create_link
from Qt.QtWidgets import QWidget, QVBoxLayout, QGridLayout, QPushButton, QLabel, QComboBox, QCheckBox, QFileDialog, QLineEdit
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import mdtraj as md
import numpy as np
import os
import ast

class AmberNetworkTool(ToolInstance):
    SESSION_ENDURING = True
    SESSION_SAVE = True

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)

        # Main container
        self.tool_window = MainToolWindow(self)

        layout = QVBoxLayout()
        entry_layout = QGridLayout()

        # File selector
        self.node_label = QLabel("Nodes Dictionary File")
        self.node_name = QLineEdit()
        self.node_name.setMinimumWidth(50)
        self.node_name.setPlaceholderText("select file")
        btn_node_file = QPushButton("Choose file (.txt)")
        btn_node_file.clicked.connect(self.choose_node_file)
        entry_layout.addWidget(self.node_label, 0, 0)
        entry_layout.addWidget(self.node_name, 0, 1)
        entry_layout.addWidget(btn_node_file, 0, 2)
        
        self.edge_label = QLabel("Edges Dictionary File")
        self.edge_name = QLineEdit()
        self.edge_name.setMinimumWidth(50)
        self.edge_name.setPlaceholderText("select file")
        btn_edge_file = QPushButton("Choose file (.txt")
        btn_edge_file.clicked.connect(self.choose_edge_file)
        entry_layout.addWidget(self.edge_label, 1, 0)
        entry_layout.addWidget(self.edge_name, 1, 1)
        entry_layout.addWidget(btn_edge_file, 1, 2)
        layout.addLayout(entry_layout)

        # File selector
#        self.out_label = QLabel("MDiGest Output Folder:")
#        self.out_name = QLineEdit()
#        self.out_name.setMinimumWidth(50)
#        self.out_name.setPlaceholderText("select MDiGest folder")
#        btn_out_file = QPushButton("Choose Folder")
#        btn_out_file.clicked.connect(self.choose_out_file)
#        entry_layout.addWidget(self.out_label, 2, 0)
#        entry_layout.addWidget(self.out_name, 2, 1)
#        entry_layout.addWidget(btn_out_file, 2, 2)
#        layout.addLayout(entry_layout)
        
#        param_layout = QGridLayout()
#        self.root1_name_label = QLabel("File Name Root #1:")
#        self.root1_name = QLineEdit()
#        self.root1_name.setPlaceholderText("required")
#        param_layout.addWidget(self.root1_name_label, 0, 0)
#        param_layout.addWidget(self.root1_name, 0, 1)
#        self.root2_name_label = QLabel("File Name Root #1")
#        self.root2_name = QLineEdit()
#        self.root2_name.setPlaceholderText("for subtraction of apo set")
#        param_layout.addWidget(self.root2_name_label, 1, 0)
#        param_layout.addWidget(self.root2_name, 1, 1)


        # Analysis type
 #       anal_layout = QGridLayout()
 #       self.combo_label = QLabel("Correlation Type:")
 #       self.combo = QComboBox()
 #       self.combo.addItems(["--Select--", "PCC", "DCC", "Covar_Disp"])
 #       param_layout.addWidget(self.combo_label, 2, 0)
 #       param_layout.addWidget(self.combo, 2, 1)
 #       layout.addLayout(param_layout)
       
        # Compute button
        btn_compute = QPushButton("Compute")
        btn_compute.clicked.connect(self.compute)
        layout.addWidget(btn_compute)

        # Attach layout to tool window
        self.tool_window.ui_area.setLayout(layout)

        self.tool_window.manage(placement="side")
        # Trajectory file placeholder
        self.directory = None

#    def choose_out_file(self):
#        directory = QFileDialog.getExistingDirectory(self, "Select MDiGest Output Directory")
#        if directory:
#            self.directory = directory
#            self.out_name.setText(directory)

    def choose_node_file(self):
        fname, _ = QFileDialog.getOpenFileName(None, "Open Nodes Dictionary", "", "Nodes Files (*.txt);;All Files (*)")
        if fname:
            self.nodefile = fname
            self.node_name.setText(fname)
 
    def choose_edge_file(self):
        fname, _ = QFileDialog.getOpenFileName(None, "Open Edges Dictionary", "", "Edges Files (*.txt);;All Files (*)")
        if fname:
            self.edgefile = fname
            self.edge_name.setText(fname)
    
#    def choose_parm_file(self):
#        fname, _ = QFileDialog.getOpenFileName(None, "Open Nodes Dictionary", "", "Edges Files (*.txt);;All Files (*)")
#        if fname:
#            self.top = fname
#            if fname.split(".")[-1] == "pdb":
#                run(self.session, f"open {fname}")
#            else:
#                import #tempfile
#                import parmed as pmd
#                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as temp_file:
#                    temp_filepath = temp_file.name
#                amber = pmd.load_file(fname, self.trajfile)
#                amber.coordinates = amber.coordinates[:,:]
#                amber.save(temp_filepath, format="pdb", charmm=True, overwrite=True, renumber=True)
#                run(self.session, f"open {temp_filepath}")
#            self.parm_name.setText(fname)

    def compute(self):
        #text = self.combo.currentText()

        #topology_0 = self.top
        #trajectory_0 = self.trajfile

        #if self.root1_name.text().strip():
        #    path1 = os.path.join(self.directory, self.root1_name.text().strip())
        #else:
        #    self.session.logger.warning("Please enter the first file root")

        #if self.root2_name.text().strip():
        #    path2 = os.path.join(self.directory, self.root2_name.text().strip())
        #else:
        #    self.session.logger.warning("You have not entered a second Filename Root, analysis will only be done on the first dataset")

        #first_load = sd.MDSdata()
        #first_load.load_from_file(file_name_root=path1)
        #if path2 is not None:
        #    second_load = sd.MDSdata()
        #    second_load.load_from_file(file_name_root=path2)

        #ss = MDS_analysis()
        #ss.load_system(topology_0, trajectory_0)
        #ss.set_selection('protein', sys_str_sel='protein')
        #ss.do_ss_calculation()

        #if text == "DCC":
        #    gcc = first_load.dcc_allreplicas['rep_0'] - second_load.dcc_allreplicas['rep_0']
        #elif text == "PCC":
        #    gcc = first_load.pcc_allreplicas['rep_0'] - second_load.pcc_allreplicas['rep_0']
        #elif text == "Covar_Disp":
        #    gcc = first_load.covar_disp_allreplicas['rep_0'] - second_load.covar_disp_allreplicas['rep_0']
        #else:
        #    self.session.logger.warning("Please select a correlation function to determine network")

        #rd = self.directory
        #ss_network(ss_stats=ss.ss_stats,
        #          gcc = gcc,
        #          nodes_save_path = rd + f'nodes_{text}.txt',
        #          edges_save_path = rd + f'edges_{text}.txt',
        #          num_sd=7)

        mdigest_network_txt(self.session, self.nodefile, self.edgefile)



def mdigest_network_txt(session, nodes_file, edges_file):
    # --- read text file ---
    with open(nodes_file, 'r') as f:
        nodes_str = f.read()
    nodes = ast.literal_eval(nodes_str)

    with open(edges_file, 'r') as f:
        edges_str = f.read()
    edges = ast.literal_eval(edges_str)

    # --- make a marker set for visualization ---
    mset = MarkerSet(session, name="Residue Network")
    session.models.add([mset])

    node_markers = {}
    # place a marker at the centroid of residues for each node
    for nid, res_list in nodes.items():
        coords = []
        run(session, f"show :{res_list} surface")
        run(session, f"color :{res_list} cornflowerblue")
        for rnum in res_list:
            try:
                model = session.models.list()[0]
                residue = [r for r in model.residues if r.number == rnum][0]
                atom = residue.find_atom("CA")
                coords.append(atom.coord)
                run(session, f"transparency #1 50")
                
            except Exception as e:
                session.logger.warning(f"{e}")
                continue
        if coords:
            centroid = tuple(np.mean(coords, axis=0))
            node_markers[nid] = mset.create_marker(centroid, (241, 0, 132, 255), 1)
            run(session, f"show #2 surface")

    # draw edges between node markers
    for key, w in edges.items():
        if isinstance(key, str):  # e.g. '9, 56'
            i, j = map(int, key.split(","))
        else:
            i, j = key
        if i in node_markers and j in node_markers:
            radius = max(0.05, abs(w)/10)
            color = (255, 0, 0, 255) if w > 0 else (0, 223, 98, 255)
            try:
                create_link(node_markers[i], node_markers[j], rgba=color, radius=radius)
            except TypeError:
                pass

    #run(session, f"color #2 blue")
