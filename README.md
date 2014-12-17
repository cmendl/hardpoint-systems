Molecular dynamics simulations of classical anharmonic hard-point particle chains
=================================================================================

C source code (with minor improvements) and Mathematica test files for the calculations in Ref. 1 and 2 below. The code uses the FFTW library and consists of three subprojects:
- monoatomic chains with a shoulder-shaped interaction potential
- biatomic hard-point particle chains with alternating masses and pairwise elastic collisions at zero separation distance
- biatomic hard-point particle chains with alternating masses and square-well interaction potential; i.e., the particles additionally "collide" when reaching a maximum separation distance

How to compile the source code:
- Windows: Visual Studio project files are provided in the *vcproj* subfolders
- Linux/Unix-type OS: makefiles are available in the *bin* and *test* subfolders. You probably have to adapt the FFTW library paths to your local installation. Some test files use a smaller lattice size (e.g., 16 instead of 4096) and require a temporary change of the 'NUM_SITES' preprocessor constant in `field.h` prior to compilation.


License
-------
Copyright (c) 2014, Christian B. Mendl  
All rights reserved.  
http://christian.mendl.net

This program is free software; you can redistribute it and/or
modify it under the terms of the Simplified BSD License
http://www.opensource.org/licenses/bsd-license.php


References
----------
1. Christian B. Mendl, Herbert Spohn  
   Current fluctuations for anharmonic chains in thermal equilibrium  
   [arXiv:1412.4609](http://arxiv.org/abs/1412.4609)
2. Christian B. Mendl, Herbert Spohn  
   Equilibrium time-correlation functions for one-dimensional hard-point systems  
   Phys. Rev. E 90, 012147 (2014), [arXiv:1403.0213](http://arxiv.org/abs/1403.0213)
3. Christian B. Mendl, Herbert Spohn  
   Dynamic correlators of Fermi-Pasta-Ulam chains and nonlinear fluctuating hydrodynamics  
   Phys. Rev. Lett. 111, 230601 (2013), [arXiv:1305.1209](http://arxiv.org/abs/1305.1209)
4. Herbert Spohn  
   Nonlinear fluctuating hydrodynamics for anharmonic chains  
   J. Stat. Phys. 154, 1191-1227 (2014), [arXiv:1305.6412](http://arxiv.org/abs/1305.6412)
