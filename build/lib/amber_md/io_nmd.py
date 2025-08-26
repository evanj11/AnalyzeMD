import numpy as np

def read_nmd(path):
    """
    Parse a .nmd file and return a dict with:
      - atom_indices
      - vectors (n_modes x n_atoms x 3)
      - frequencies
    """
    atom_indices = []
    frequencies = []
    vectors = []

    with open(path, "r") as f:
        lines = f.readlines()

    mode_data = []
    for line in lines:
        parts = line.strip().split()
        if not parts:
            continue
        tag = parts[0].lower()
        if tag == "mode":
            # Save previous mode
            if mode_data:
                vectors.append(np.array(mode_data))
                mode_data = []
            freq = float(parts[2]) if len(parts) >= 3 else None
            frequencies.append(freq)
        elif tag == "resids":
            # Atom line: atom index + 3 components
            idx = int(parts[1])
            disp = [float(x) for x in parts[2:5]]
            atom_indices.append(idx)
            mode_data.append(disp)

        print(frequencies)
    # Last mode
    if mode_data:
        vectors.append(np.array(mode_data))

    vectors = np.array(vectors)  # shape: (n_modes, n_atoms, 3)

    return {
        "atom_indices": np.unique(atom_indices).tolist(),
        "vectors": vectors,
        "frequencies": frequencies,
    }

