from chimerax.core.errors import UserError
import mdtraj as md
import numpy as np

def open_nc(session, path, file_name, **kw):
    try:
        traj = md.load(path)
    except Exception as e:
        raise UserError(f"Could not read NetCDF trajectory: {e}")

    from chimerax.atomic import AtomicStructure
    model = AtomicStructure(session, name=file_name)
    session.logger.info(f"Loaded trajectory with {traj.n_frames} frames:, {traj.n_atoms} atoms")
    return [model], None

def save_nc(session, path, models, **kw):
    import mdtraj as md
    from mdtraj.formats import NetCDFTrajectoryFile

    model = models[0]
    coords = model.positions.array
    with NetCDFTrajectoryFile(path, 'w', force_overwrite=True) as f:
        f.write(coords[np.newaxis])
    session.logger.info(f"Saved trajectory to {path}")

