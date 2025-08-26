# Minimal NMWiz (.nmd) loader + visualizer for ChimeraX
# - Parses .nmd (ProDy/NMWiz)
# - Maps site vectors (usually Cα) to model atoms
# - Shows arrows on sites, or expands to all atoms per residue for animation

import numpy as np
from chimerax.core.commands import run
from chimerax.core.models import Model
from chimerax.geometry import Place
from chimerax.ui.gui import MainToolWindow
from chimerax.core.tools import ToolInstance
from chimerax.atomic import AtomicStructure, selected_atoms
from chimerax.core.commands import run 
from Qt.QtWidgets import QWidget, QVBoxLayout, QGridLayout, QPushButton, QLabel, QComboBox, QCheckBox, QFileDialog, QLineEdit
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

# ---------- Parsing ----------

def parse_nmd(path):
    """Parse an .nmd file. Returns dict with fields:
       coordinates (N x 3) optional
       atomnames, resnames, chainids (list[str], length N)
       resids (list[int], length N)
       modes: dict[int -> (N x 3) ndarray]
    """
    data = {"modes": {}}
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            tag = parts[0].lower()
            vals = parts[1:]

            if tag == "coordinates":
                arr = np.array([float(x) for x in vals], dtype=float)
                data["coordinates"] = arr.reshape(-1, 3)
            elif tag == "atomnames":
                data["atomnames"] = vals
            elif tag == "resnames":
                data["resnames"] = vals
            elif tag == "chainids":
                data["chainids"] = vals
            elif tag == "resids":
                # resids can be ints; keep as strings *and* ints for robust matching
                try:
                    data["resids"] = [int(x) for x in vals]
                except ValueError:
                    data["resids"] = vals
            elif tag == "mode":
                mode_index = int(vals[0])
                arr = np.array([float(x) for x in vals[2:]], dtype=float)
                vec = arr.reshape(-1, 3)
                data["modes"][mode_index] = vec

    # Basic sanity
    # Determine Nsites from whichever array we have
    n_candidates = []
    for key in ("coordinates", "atomnames", "resnames", "chainids", "resids"):
        if key in data:
            if key == "coordinates":
                n_candidates.append(len(data[key]))
            else:
                n_candidates.append(len(data[key]))
    if n_candidates:
        N = n_candidates[0]
        for m, vec in data["modes"].items():
            if len(vec) != N:
                raise ValueError(f".nmd mode {m} length {len(vec)} != Nsites {N}")
    return data

# ---------- Mapping to model ----------

def build_site_keys(nmd, selection):
    """Keys for matching: (chainID, resid, atomname). If atomnames not present, use 'CA'."""
    chainids = nmd.get("chainids")
    resids   = nmd.get("resids")
    atomnames = nmd.get("atomnames", None)
    if chainids is None or resids is None:
        raise ValueError("NMD missing chainids/resids; cannot map.")
    # Normalize resids to str for comparison, but keep ints to try both
    keys = []
    for i in selection:
        if len(chainids) == len(resids):
            chain = str(chainids[i])
        else:
            chain = str("A")
        resid = resids[i]
        resid_str = str(resid)
        atom = atomnames[i] if atomnames else "CA"  # default to CA
        keys.append( (chain, resid_str, atom) )
    return keys

def model_atom_index_map(model):
    """Return dict keyed by (chainID, resid_str, atomname) -> atom object."""
    mp = {}
    for a in model.atoms:
        r = a.residue
        chain = getattr(r, "chain_id", "") or ""     # '' for no chain
        resid_str = str(getattr(r, "number", r.number))
        key = (chain, resid_str, a.name)
        mp[key] = a
    return mp

def model_residue_atom_lists(model):
    """Map (chainID, resid_str) -> list of atoms in that residue."""
    res_map = {}
    for r in model.residues:
        chain = getattr(r, "chain_id", "") or ""
        resid_str = str(getattr(r, "number", r.number))
        res_map[(chain, resid_str)] = list(r.atoms)
    return res_map

# ---------- Visualization ----------

def show_site_arrows(session, model, nmd, selection, color="red", mode_index=1, scale=3.0):
    # Create a container model for arrows
    arrow_group = Model("Mode Arrows", session)
    session.models.add([arrow_group], parent=model)

    coords = model.atoms.coords

    # compute molecular size for scaling
    minc, maxc = coords.min(axis=0), coords.max(axis=0)
    mol_size = ((maxc - minc)**2).sum()**0.5
    arrow_scale = scale * mol_size

    mode = nmd["modes"][mode_index]
    keys = build_site_keys(nmd, selection)
    amap = model_atom_index_map(model)
    
    matched = 0
    for vec, key in zip(mode, keys):
        atom = amap.get(key)
        if atom is None:
            continue  # site not present in model
        start = atom.coord
        end = start + arrow_scale * vec

        run(session, f"shape cylinder fromPoint {start[0]},{start[1]},{start[2]} toPoint {end[0]},{end[1]},{end[2]} radius 0.2 color {color}", log=False)

        # Grab the most recently created model (the cylinder)
        cyl = session.models.list()[-1]
        # Re-parent it under the arrow group
        session.models.remove([cyl])
        arrow_group.add([cyl])

        cone_end = end + 20*vec
        run(session, f"shape cone fromPoint {end[0]},{end[1]},{end[2]} toPoint {cone_end[0]},{cone_end[1]},{cone_end[2]} radius 0.3 color {color}", log=False)
        cone = session.models.list()[-1]
        # Re-parent it under the arrow group
        session.models.remove([cone])
        arrow_group.add([cone])

    return arrow_group

def expand_mode_to_all_atoms(model, nmd, mode_index=1):
    """Expand site vectors (e.g., Cα) to all atoms in each residue (translation only).
       Returns an (N_atoms x 3) array aligned to model.atoms order.
    """
    mode = nmd["modes"][mode_index]     # shape (Nsites, 3)
    keys = build_site_keys(nmd)         # (chain, resid_str, atomname)
    res_lists = model_residue_atom_lists(model)

    # Build residue-level vector: use the site's (chain,resid) vector; if multiple sites per residue,
    # we average (rare for CA-only).
    from collections import defaultdict
    res_vec_sum = defaultdict(lambda: np.zeros(3, dtype=float))
    res_vec_cnt = defaultdict(int)

    for vec, (chain, resid_str, atomname) in zip(mode, keys):
        res_key = (chain, resid_str)
        res_vec_sum[res_key] += vec
        res_vec_cnt[res_key] += 1

    res_vec = {k: (res_vec_sum[k] / res_vec_cnt[k]) for k in res_vec_sum}

    # Now make per-atom array
    all_atoms = list(model.atoms)
    out = np.zeros((len(all_atoms), 3), dtype=float)
    for idx, atom in enumerate(all_atoms):
        r = atom.residue
        chain = getattr(r, "chain_id", "") or ""
        resid_str = str(getattr(r, "number", r.number))
        vec = res_vec.get((chain, resid_str))
        if vec is not None:
            out[idx] = vec  # translation for all atoms in residue
        # else remains zero
    return out

# ---------- Convenience wrappers ----------

def nmd_show_sites(session, mode_path, selection, color="red", mode_index=1, scale=3.0):
    """Public: show arrows on NMD sites (e.g., CA)."""
    nmd = parse_nmd(mode_path)
    # use first AtomicStructure model
    ms = [m for m in session.models if m.__class__.__name__ == "AtomicStructure"]
    if not ms:
        session.logger.warning("No atomic model open.")
        return
    model = ms[0]
    show_site_arrows(session, model, nmd, selection, color, mode_index=mode_index, scale=scale)


class AmberNMWTool(ToolInstance):
    SESSION_ENDURING = True
    SESSION_SAVE = True

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)

        # Main container
        self.tool_window = MainToolWindow(self)

        layout = QVBoxLayout()
        entry_layout = QGridLayout()

        # File selector
        self.nmd_label = QLabel("Normal Mode File")
        self.nmd_name = QLineEdit()
        self.nmd_name.setMinimumWidth(50)
        self.nmd_name.setPlaceholderText("select file")
        btn_nmd_file = QPushButton("Choose NMD (.nmd)")
        btn_nmd_file.clicked.connect(self.choose_nmd_file)
        entry_layout.addWidget(self.nmd_label, 0, 0)
        entry_layout.addWidget(self.nmd_name, 0, 1)
        entry_layout.addWidget(btn_nmd_file, 0, 2)

        self.sel_label = QLabel("Selection")
        self.sel_entry = QLineEdit()
        self.sel_entry.setPlaceholderText("select on model atoms for NMD")
        entry_layout.addWidget(self.sel_label, 1, 0)
        entry_layout.addWidget(self.sel_entry, 1, 1)

        self.color_label = QLabel("Color of Vector Arrows")
        self.color_text = QLineEdit()
        self.color_text.setPlaceholderText("default is red")
        entry_layout.addWidget(self.color_label, 2, 0)
        entry_layout.addWidget(self.color_text, 2, 1)

        layout.addLayout(entry_layout)

        btn_compute = QPushButton("Compute")
        btn_compute.clicked.connect(self.compute)
        layout.addWidget(btn_compute)

        self.tool_window.ui_area.setLayout(layout)
        self.tool_window.manage(placement="side")

    def choose_nmd_file(self):
        fname, _ = QFileDialog.getOpenFileName(None, "Open Normal Mode Data", "", "NMD File (*.nmd);;All Files (*)")
        if fname:
            self.nmdfile = fname
            self.nmd_name.setText(fname)

    def compute(self):
        if not self.nmdfile:
            self.session.logger.warning("No NMD file selected.")
            return
                
        chains = self.session.models.list()[0].chains
        if len(chains) <= 1:
            run(self.session, "changechains #1 A")
            

        sel = selected_atoms(self.session)
        if len(sel) > 0:
            selection = [a.coord_index for a in sel]
            self.sel_entry.setText("Using selected")
        else:
            self.session.logger.warning("Please select atoms to display normal mode data")

        nmd_show_sites(self.session, self.nmdfile, selection, color=self.color_text.text().strip(), mode_index=1, scale=3.0)

