function y=osgood1(x)
global Kaxial Naxial Eyoung sapoisson
y=x/Eyoung+(x/Kaxial)^(1/Naxial)-sapoisson;