%This script finds a random trial point and writes the result into the 
%'WorkingDirectory'. Also updates the Array of points tried out

%!! This needs to be updated so that it uses a metamodel and optimizes it at
%every iteration !!
function [ChangedODMatrix,rowIDs,currTextFilename] = FindTrialPoint(iter,baseODMatrix,AllowedReductionPercentage,numNewStations,HOMEDIRECTORY)

N = size(baseODMatrix,1);

rowIDs = randperm(N);
rowIDs = rowIDs(1:numNewStations);

ChangedODMatrix = baseODMatrix;
ChangedODMatrix(rowIDs,3)=ChangedODMatrix(rowIDs,3)*(1-AllowedReductionPercentage);
currTextFilename = [HOMEDIRECTORY '\\TrialPoints\\Iter_' num2str(iter) '.txt'];
dlmwrite(currTextFilename,ChangedODMatrix,'\t');
end

