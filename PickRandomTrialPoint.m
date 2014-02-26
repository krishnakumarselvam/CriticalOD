%This function picks and returns a random point from the file
%'Inputs/RandomTrialPoints.mat 

%A trial point is defined by the (assumed) ordered list of OD Demand
%values.

function [CurrTrialPoint] = PickRandomTrialPoint()
load ('Inputs/RandomTrialPoints', 'randomTrialPoints');
n = randi(size(randomTrialPoints,1));
CurrTrialPoint = randomTrialPoints(n,:);
end