function [xStar,f_sim_xk_TT] = get_fSim1_TT(fileX,fileNwk,fileQgParam,fileOut_txt,fileIn_py,file_gtOpt,currPtId,allObsFile,allPtsFile,runAimsun,loadGtsFromQg_bool,currGT_notFromQg,fileAllRoundedGTs,dbnIter,initPtId,dirName)
% differs from get_fSim.m : dirName added to output files
%
% fileX:        contains x_out to load. ie x*
% fileNwk:      nwk params: subNwk11cData*.mat
% fileOut_txt:  contains the old and new signal plans (before and after rounding)
% fileIn_py:    python file with the controlDico variable (new set of GTs)
% file_gtOpt:   mat file with the gtOpts corresponding to controlDico
% currPtId:   ID of the current iteration (ie k, and not k-1): this will save CentreV5 with the new GT plan as the model: CentreV5_k.ang that will be used to calc f_sim for this iteration.
% adapt plan (e.g. rounding may be necessary)
% write python script to integrate it within the model

if runAimsun ~= 0
    
    if loadGtsFromQg_bool
        load(fileX,'x_out');
        load(fileQgParam,'last_index');
        load(fileNwk,'cycleTimes_s','phase2Node');
        nbPhases = length(phase2Node);
        xStar = x_out(last_index:(last_index+nbPhases-1));
        clear x_out last_index cycleTimes_s phase2Node;
    else
        xStar = currGT_notFromQg;
    end;
    
    gt2python10(fileX,fileNwk,fileQgParam,fileOut_txt,fileIn_py,file_gtOpt,currPtId,loadGtsFromQg_bool,currGT_notFromQg,dbnIter,initPtId,dirName);
    % may differ from x_out because x_kUsed may be rounded
    
    
    newAvgTTFile = angParams2python10(currPtId,dbnIter,initPtId,dirName);
    
    
    %{
 writes:
- the new model name
- the number of replications
- the new avgTT observation
in a python file so that it can be imported when calling other python files (e.g. integrateNewGT_mdl.py)
passing them as arguments was too difficult via matlab
- returns a string containing the name of the new file that contains the avgTT's of each replication of those runs
    %}
    
    cd 'C:/Program Files/TSS-Transport Simulation Systems/AIMSUN 6.1/'
    % from the 'AIMSUN NG v5.1' directory we can then call:
    % we need to be in that directory or else the python library files are not recognized by angconsole.exe (pb with the new vsn of aimsun but I stopped trying to debug it)
    % NOTE: current working directory has changed
    
    !aconsole.exe -script C:/MIT/Lab/KanchanaNanduri/CodeSO_TT_SD/Runs1/PyInputs/integrateNewGT_mdl_newReplic.py
    % this opens the model tmp.ang
    % prints out the new and the old plans: TR/diffPlans_aimsunOutput.txt
    % creates a new model TR/Aimsun/CentreV5_currPtId.ang with the new
    % plans.
    
    disp('launching aimsun run');
    
    while ~exist(newAvgTTFile)
        pause(10) % pauses for 10 sec.
        !aconsole.exe -script C:/MIT/Lab/KanchanaNanduri/CodeSO_TT_SD/Runs1/PyInputs/simRepAll_newMdl_newReplic.py 
        %>> out1_115d.txt
    end;
    % runs the new mdl, with nbReplic replications, stores the avgTT of each replication in newAvgTTFile
    % after each replication the mdl is no longer saved (because the same model is run with different seeds, and for now there is no need to save each model)
    % checked that the avgTT calculations are coherent with the data in vehData
    % and that in vehData only exits after the warmUpPeriod were considered (ie <= 900sec ie 15min)
    
    disp('aimsun run finished');
    
    cd 'C:/MIT/Lab/KanchanaNanduri/CodeSO_TT_SD/Runs1/';
    
    replicID = getReplicID10(dirName);
    
else
    replicID = -1;
    newAvgTTFile = rand(1,1);
    xStar = get_rdmGT_uniformly(fileNwk);
    file_gtOpt = round(xStar);
end;



addRoundedGTs(file_gtOpt,fileAllRoundedGTs,replicID,currPtId);
f_sim_xk_TT = addObservtns(newAvgTTFile,currPtId,allObsFile,runAimsun,replicID);
% one file per plan: read it into matlab, get the avg per replications,
% write that into file allObsFile: each line := pointID (must be unique), avgDelay

addPoint(xStar,currPtId,allPtsFile,replicID);
% write the current point (ie decision vble) into file allPointsFile:
% each line := pointID, vector of values (eg x1, ..., xDim; where Dim is the dimension of the point, eg the number of endogenous phases)
