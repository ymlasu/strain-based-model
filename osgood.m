function y=osgood(x)
global Kaxial Naxial Eyoung sa
y=x/Eyoung+(x/Kaxial)^(1/Naxial)-sa;