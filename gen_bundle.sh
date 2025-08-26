#!/bin/bash

rm -rf build dist ChimeraX_AmberMD.egg-info   

${CHIMERAX_EXE} --nogui --cmd "devel install . ; exit"
