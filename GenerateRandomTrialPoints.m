%This script generates a random set of trial points from which we can pick
%one to start the algorithm. It's useful in two ways. One is that we don't
%have to spend the computational effort in generating random trial points
%every time we run the algorithm. The other is that we can for instance
%start with the same trial point every time, and run controlled
%experiments.

%Logic:
%We need to remove a certain number of vehicles from the network. We do this by
%randomly choosing a few OD pairs, and then remove 25 % of the demand from
%each of the OD pair.

%In the end this script just saves the third column of the OD pair table.
%It assumes that the order of OD pairs is maintained in the rest of the
%script.

%Parameters to set: 

%From now on, only the Demand values belonging to the TopODIndices OD pairs
%will be stored in the variables throughout the algorithm. However, the
%text file that is written for Aimsun Model creation willcontinue to have
%all the OD pair values.
clc
clear
FRACTION_TO_REDUCE    = 0.25;
NUM_VEHICLES_TO_REMOVE  = 40; %Change this in the main file as well
NUM_PTS_TO_GENERATE     = 100; % We generate and store these many random trial points 

baseODMatrix = textread('Inputs/ODpairs.txt');%Read from the Original OD pair data
N = size(baseODMatrix,1); %Number of OD pairs in the network
PROBLEMDIMENSION = 10;
randomTrialPoints = zeros (NUM_PTS_TO_GENERATE, PROBLEMDIMENSION); % This is where we store our results


[temp,ind] = sort(baseODMatrix(:,3),'descend');
TopODIndices = ind(1:PROBLEMDIMENSION); % This now contains the index of the OD pairs with the 200 (i.e. numODsConsidered) highest demands.

for i = 1:NUM_PTS_TO_GENERATE
    disp(['Iteration ' num2str(i)]);
    rowIDs = randperm(PROBLEMDIMENSION);
    ChangedODMatrixList = baseODMatrix(:,3)';
    vehiclesRemoved = 0;
    
    j=1;
    while (vehiclesRemoved <NUM_VEHICLES_TO_REMOVE)
        temp = vehiclesRemoved + ChangedODMatrixList(TopODIndices(rowIDs(j)))*(1-FRACTION_TO_REDUCE);
        if(temp>=NUM_VEHICLES_TO_REMOVE)
            currentVehRemoved = NUM_VEHICLES_TO_REMOVE - (vehiclesRemoved);
        else
            currentVehRemoved = ChangedODMatrixList(TopODIndices(rowIDs(j)))*(1-FRACTION_TO_REDUCE);
        end
        vehiclesRemoved = vehiclesRemoved + currentVehRemoved;
        ChangedODMatrixList(TopODIndices(rowIDs(j))) = ChangedODMatrixList(TopODIndices(rowIDs(j))) - currentVehRemoved;
        j=j+1;
    end
    randomTrialPoints(i,1:PROBLEMDIMENSION) = ChangedODMatrixList(TopODIndices);
end

save ('Inputs/RandomTrialPoints', 'randomTrialPoints','TopODIndices');