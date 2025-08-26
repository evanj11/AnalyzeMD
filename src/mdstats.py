from chimerax.core.commands import CmdDesc, register, StringArg, BoolArg
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def compute_mdstats(session, trajfile, topfile, stat="rmsd", free_energy=False):
    traj = md.load(trajfile, top=topfile)
    list_stat = ["rmsd", "rmsf", "radg"]

    if stat == "rmsd":
        values = md.rmsd(traj, traj, 0)
        ylabel = "RMSD (nm)"
    elif stat == "rmsf":
        mean_coords = traj.xyz.mean(axis=0)
        values = np.sqrt(((traj.xyz - mean_coords) ** 2).sum(axis=2).mean(axis=0))
        ylabel = "RMSF (nm)"
    elif stat == "radg":
        values = md.compute_rg(traj)
        ylabel = "Radius of Gyration (nm)"
    else:
        raise ValueError(f"Unknown stat: {stat}, please choose from {list_stat}")

    session.logger.info(f"{stat.upper()} calculated: {len(values)} values")

    plt.figure()
    if stat in ["rmsd", "radg"]:
        plt.plot(values)
        plt.xlabel("Atom Index")
        plt.ylabel(ylabel)
        plt.title(f"{stat.upper()} over time")
    elif stat == "rmsf":
        plt.bar(np.arange(len(values)), values)
        plt.xlabel("Atom Index")
        plt.ylabel(ylabel)
        plt.title("RMSF per atom")

    plt.tight_layout()
    plt.show(block=False)

    if free_energy:
        rmsd = md.rmsd(traj, traj, 0)
        rg = md.compute_rg(traj)
        
        H, xedges, yedges = np.histogram2d(rmsd, rg, bins=50)
        F = -np.log(H + 1e-12)
        F -= F.min()

        plt.figure()
        plt.contourf(xedges[:-1], yedges[:-1], F.T, levels=50, cmap='seismic')
        plt.xlabel("RMSD (nm)")
        plt.ylabel("Rg (nm)")
        plt.title("Free Energy Landscape")
        cbar = plt.colorbar()
        cbar.set_label("ΔG (kBT)")
        plt.tight_layout()
        plt.show(block=False)

    return values

desc = CmdDesc(
    required=[("trajfile", StringArg), ("topfile", StringArg)],
    optional=[("stat", StringArg)],
    keyword=[("free_energy", BoolArg)],
)

#def register_command(logger):
#    register("mdstats", desc, compute_mdstats, logger=logger)
