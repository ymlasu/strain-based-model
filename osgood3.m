function y=osgood3(x)
global Kaxial Naxial Eyoung Strain
y=x/Eyoung+(x/Kaxial)^(1/Naxial)-Strain;