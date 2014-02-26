%This script finds a new trial point and writes the result into the 
%'WorkingDirectory'. Also updates the Array of points tried out

%In this version, the script optimizes a metamodel and finds the next trial
%point. Also, I use the updated formulation, in which we can remove a fixed
%number of carrs from the network, and we need to pick how many from each
%of the OD pairs.

function [TrialPoint] = FindTrialPoint(iter,baseODMatrix,HOMEDIRECTORY)

if(iter <= 100)
    %We need to select a random trial point;
    TrialPoint = PickRandomTrialPoint;
    
else
    %We have to optimize the current metamodel to get the new trial point.
end
ChangedODMatrix = baseODMatrix;
ChangedODMatrix(:,3) = TrialPoint';
currTextFilename = [HOMEDIRECTORY '\\TrialPoints\\Iter_' num2str(iter) '.txt'];
dlmwrite(currTextFilename,ChangedODMatrix,'\t');

end

