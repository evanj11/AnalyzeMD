CHIMERAX_BUILD=$(which chimerax_build)

# build and install from the correct folder
$CHIMERAX_BUILD build src/amber_md
$CHIMERAX_BUILD install src/amber_md
$CHIMERAX_BUILD dist src/amber_md
