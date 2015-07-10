RivetMC_BOOSTEDBOSON.so: MC_BOOSTEDBOSON.cc
	rivet-buildplugin RivetMC_BOOSTEDBOSON.so MC_BOOSTEDBOSON.cc `fastjet-config --prefix`/lib/libEnergyCorrelator.a
