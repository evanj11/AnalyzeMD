# AnalyzeMD
ChimeraX plug-in for analysis of molecular dynamics simulations

## Install
To install the AnalyzeMD plug-in for ChimeraX, you first must install the dependencies:
- For MacOS/Linux CLI
  
  ```
  export CHIMERAX_EXE=/path/to/chimerax
  ${CHIMERAX_EXE} -m pip install mdtraj deeptime scikit-learn parmed
  bash gen_bundle.sh
  ```
 
- For Windows PowerShell
  
  ```
  & "C:\Program Files\ChimeraX\bin\ChimeraX.exe" -m pip install mdtraj deeptime scikit-learn parmed
  ```
  
  and after successful installation of dependencies, within ChimeraX:

  ```
  devel install PATH_TO_SOURCE_CODE_FOLDER
  ```

## Usage
AnalyzeMD is used to automate several beginning-stages analysis steps for Molecular Dynamics Simulations and includes four modules:
1. **Simulation Characterization**
  - Allows for the calculation of RMSD/RMSF/RadG across the simulation on the whole protein or a specified selection\textsuperscript{\ast}
  - Can generate the Free Energy Landscape (FEL) of a given selection with an option to input the temperature for $\Delta G$ determination
  - Can calculate pairwise graphs such as RMSD, Principal Component Analysis (first two PCs), and Dynamic Cross Correlation Matrices
2. **Geometry Analysis**
  - Plots the distance between two selected atoms or the angle between three selected atoms across the simulation
  - Determines the number of H-Bonds within a given selection, including a contact map option, throughout the simulation
3. **Normal Mode Wizard**
  - Graphs the first normal mode eigenvectors from a .nmd file onto the displayed model
  - Can be limited to the currently-selected atoms
4. **Markov State Models**
  - Utilizes either RMSD/RadG or Phi-Psi couples to determine discreet kinetic states
  - User can input the number of microstates to perform KMeans clustering on\textsuperscript{\dagger} and the number of macrostates that PCCA+ coarse-grains the microstates into\textsuperscript{\ddagger}
  - The reactive flux path is then drawn over the PCCA+ macrostates to illustrate a potential path traveled

$^\ast$ All `AnalyzeMD` functions support selecting specific atoms on the model and using that as the selection<br />
$^\dagger$ As a rule of thumb, the number of microstates should be ~1/1,000 of the number of frames<br />
$^\ddagger$ If you get errors citing `Index out of range`, decrease the number of microstates

- For **Simulation Characterization** and **Geometry Analysis**, the script will automatically open the topology file, and so these analyses should be accessed prior
  to opening a model.
- `AnalyzeMD` currently supports trajectory files in .dcd, .nc, .xtc and topology file in .pdb or .parm7
  - Regardless of the topology file loaded, a pdb file will be opened in the ChimeraX session
  - *Where possible, use the .pdb as the topology to prevent ChimeraX freezing and atom index mismatches*

## Potential Caveats
- The number of frames used in the simulation must be tailored to the analysis type. For example, using >1,000 frames for **Simulation Characterization** or **Geometry Analysis** 
  is not recommended; however, >50,000 frames might be necessary to achieve the statistical sampling needed for **Markov State Models**.
- Unless absolutely necessary, waters and ions should be stripped prior to use to accelerate the analysis (specific waters can be preserved)
