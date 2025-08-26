from chimerax.core.commands import CmdDesc, register, StringArg
from chimerax.core.tools import get_singleton

# Import our tool class
from .tool import AmberMDTool


def mdstats(session, traj=None, gui=False):
    """Command entry point for AmberMD analysis."""
    if gui:
        # Open the dockable GUI (singleton ensures only one copy)
        tool = get_singleton(session, AmberMDTool, "AmberMD Tool")
        tool.display(True)  # show the panel
        session.logger.info("AmberMD GUI opened.")
    else:
        # Stub for trajectory analysis
        session.logger.info(f"Running mdstats on {traj}")


# Define the command descriptor
mdstats_desc = CmdDesc(
    optional=[('traj', StringArg)],
    keyword=[('gui', bool)],
    synopsis="Compute RMSD/RMSF/Rg or open the AmberMD GUI."
)


def register_command(logger):
    register('mdstats', mdstats_desc, mdstats, logger=logger)

