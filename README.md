# tov_solver - Solving the Tolman–Oppenheimer–Volkoff equation using the Runge-Kutta method

## Background

<img align="right" src="https://github.com/cjdietl/tov_solver/blob/master/MR.png">

This is part of the code for my paper ['Properties of gravitationally bound dark compact ultra dense objects'](http://linkinghub.elsevier.com/retrieve/pii/S0370269312001463). It basically solves the differential equations laid out in the seminal paper ['On Massive Neutron Cores'](http://prola.aps.org/abstract/PR/v55/i4/p374_1) by J. R. Oppenheimer and G. M. Volkoff. 

The Tolman-Oppenheimer-Volkoff equation allows to determine the stable regions for a neutron star - similar to the Chandrasekhar mass for white dwarfs. In our paper though, we use it with a twist: The neutron mass is replaced by the mass for a hypothetical dark matter particle, which is way higher than 'normal' baryonic matter. This results in crazy small, but very heavy objects.

The Mass-Radius relationship has a funny curl at the end (Fig.). This is actually the unstable region, where the object collapses under its own gravitational pressure and forms a black hole.

See also https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_limit .

## Running oppenheimer.cc
g++ oppenheimer.cc -o oppenheimer
./oppenheimer
