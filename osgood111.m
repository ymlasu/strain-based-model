function y=osgood111(x)
global Ktor Ntor Gshear Strain
y=x/Gshear+(x/Ktor)^(1/Ntor)-Strain;