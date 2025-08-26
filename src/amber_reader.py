import mdtraj as md
from chimerax.core.models import AtomicStructure
from chimerax.atomic import AtomicStructure as CXAtomicStructure

def open_amber(session, paths, **kw):
    """
    Open Amber parm7 + nc trajectory into ChimeraX.
    paths: list of files (parm7 + nc)
    """
    # Find parm7 and nc
    parm7 = None
    nc = None
    for p in paths:
        if p.endswith((".parm7", ".prmtop")):
            parm7 = p
        elif p.endswith(".nc"):
            nc = p
    if parm7 is None or nc is None:
        raise ValueError("Need both a .parm7 topology and .nc trajectory")

    # Load with MDTraj
    traj = md.load(nc, top=parm7)

    # Create ChimeraX AtomicStructure for first frame
    m = CXAtomicStructure(session, name="AmberMD")
    top = traj.topology

    atom_map = []
    for atom in top.atoms:
        cx_atom = m.new_atom(atom.name, atom.element.symbol if atom.element else "C")
        m.add_atom(cx_atom, atom.residue.index)
        atom_map.append(cx_atom)

    # Set coordinates for first frame
    m.set_coord(traj.xyz[0] * 10.0)  # convert nm → Å

    # Add to session
    session.models.add([m])

    # Animate with coordinate sets
    for i in range(1, traj.n_frames):
        m.add_coord_set(traj.xyz[i] * 10.0)

    return [m], f"Amber trajectory: {parm7} + {nc}"


if format_name == "Amber Parm7":
    fmt = "parm7"
elif format_name == "Amber NetCDF":
    fmt = "nc" 

fr = FileReader(
            (file_name, fmt, stream),
            just_geom=False,
            get_all=True,
            max_length=max_length,
            log=session.logger,
        )

if fr.file_type == "parm7":
    
