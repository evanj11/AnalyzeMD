# src/amber_md/__init__.py

from chimerax.core.toolshed import BundleAPI
from chimerax.core.commands import BoolArg, IntArg, ModelsArg, StringArg, register, OpenFileNameArg

class _AmberMD_API(BundleAPI):

    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        from .mdstats import mdstats
        from .mdcharact import mdcharact
        from .io_nc import io_nc
        from .amber_reader import open_amber

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        """
        register commands
        """
        if command_info.name == "mdstats":
            from .mdstats import compute_mdstats, desc
            register("mdstats", desc, compute_mdstats)

    @staticmethod
    def start_tool(session, bi, ti):
        """
        start tools
        """
        if ti.name == "Simulation Characterization":
            from .mdcharact import AmberMDTool
            tool = AmberMDTool(session, ti.name)
            return tool
        elif ti.name == "Geometry Analysis":
            from .mddist import AmberDistTool
            tool = AmberDistTool(session, ti.name)
            return tool
        elif ti.name == "Normal Mode Wizard":
            from .mdnmwiz import AmberNMWTool
            tool = AmberNMWTool(session, ti.name)
            return tool
        elif ti.name == "Markov State Model":
            from .mdmsm import AmberMSMTool
            tool = AmberMSMTool(session, ti.name)
            return tool
        elif ti.name == "Network Analysis":
            from .mdnetwork import AmberNetworkTool
            tool = AmberNetworkTool(session, ti.name)
            return tool


bundle_api = _AmberMD_API()
