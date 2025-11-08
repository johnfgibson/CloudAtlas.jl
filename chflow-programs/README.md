# CloudAtlas.jl : channelflow programs
Several channelflow programs useful for coordinating direct numerical simulations in channelflow with ODE models of CloudAtlas.jl

   1. `projectfield.cpp` Project a channelflow DNS field (`u.nc` file) onto a symmetric basis set.
   2. `projectseries.cpp` Project a time series of DNS fields onto a symmetric basis set.

These programs can be compiled and linked to the channelflow-2.0 libraries using the included Makefile. 