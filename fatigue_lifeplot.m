
% Predicted Fatigue life data
yyy = ([4.05516710226898;4.05516710226898;
    5.08154346360181;5.08764762353692;4.00718226046569;
    4.00264138941051;5.03373498432192;5.03943113129950]);

%Experimental Fatigue life data
xxx = [3.7084209;
3.792391689;
4.950909815;
5;
4.04453976;
3.991226076;
5.301029996;
5.311753861];

plot(xxx,yyy,'o')
title('304 box')
xlabel('Experimental Fatigue life logN')
ylabel('Predicted Fatigue life logN')
grid
