# AnalyzeMD
ChimeraX plug-in for analysis of molecular dynamics simulations

## Install
To install the AnalyzeMD plug-in for ChimeraX, you first must install the dependencies:
  - For MacOS/Linux
  
  ```
  export CHIMERAX_EXE=/path/to/chimerax
  ${CHIMERAX_EXE} -m pip install mdtraj deeptime scikit-learn
  bash gen_install.sh
  ```
  
  - For Windows PowerShell
  
  ```
  & "C:\Program Files\ChimeraX\bin\ChimeraX.exe" -m pip install mdtraj deeptime scikit-learn
  ```
  
  and after successfult installation of dependencies, within ChimeraX:

  ```
  devel install PATH_TO_SOURCE_CODE_FOLDER
  ```

## Usage
