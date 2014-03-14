function [TrialPoint] = FindMixturePoint (iter,baseODMatrix,HOMEDIRECTORY,TopODIndices)
    TrialPoint = PickRandomTrialPoint();
    ChangedODMatrix = baseODMatrix;
    ChangedODMatrix(TopODIndices,3) = TrialPoint';
    currTextFilename = [HOMEDIRECTORY '\\TrialPoints\\Iter_' num2str(iter) '.txt'];
    dlmwrite(currTextFilename,ChangedODMatrix,'\t');
    
    
end