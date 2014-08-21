function gt2python10(fileData,fileNwk,fileQgParam,fileOut_txt,fileOut_py,gtOptFile,currPtId,loadGtsFromQg_bool,currGT_notFromQg,dbnIter,initPtId,dirName)

% writes vector xGT of green times to a python file, that will be used to integrate these times into the AIMSUN model

% input:  matlab new traffic signal settings
%           data file gtOptFile: contains the values of the final GTs (may differ from optimal GTs due to rounding)
% output :  a python script with new traffic signal settings


%%%%%
% Processing the matlab solution
%%%%%

%{
% expl:
 fileNwk = 'GTplans/subNwk11cData_1_0.mat';
 fileData = 'GTplans/data_totTime4CVGD_fmincon_subNwk11c_muHatLargerLin_totalTimeObjFctn_piFsolve_mutFsolve_e007.mat';
 fileOut_txt = 'GTplans/tmp_newJunction_beforeAfterRounding.txt';
 fileOut_py = 'GTplans/tmp_new_gt.txt';

 load(fileNwk,'initialGreenTimePerPhase_s','phaseInfo','signalizedNode_s','fixedPhases_T_s','totNbPhasesPerNode_s','interRedGreenAfterPh_s','availableGreenTimes_s','phase2Node','cycleTimes_s');
 load(fileData,'x_out');

% for nwk11c:
% if only signalized mu's are endog: last_index = 2196;
% if all mu is endog: last_index = 2238;
% nbPhases = length(phase2Node);
% x_out(last_index:(last_index+nbPhases-1)).*cycleTimes_s(phase2Node)
%}

beforeRounding=1;
if loadGtsFromQg_bool
    print_gt_txt(fileOut_txt,fileData,fileNwk,fileQgParam,beforeRounding,[]);
else
    print_gt_txt_notFromQg(fileOut_txt,fileNwk,beforeRounding,currGT_notFromQg);
end;


load(fileNwk,'initialGreenTimePerPhase_s','signalizedNode_s','phaseInfo','availableGreenTimes_s','fixedPhases_T_s',...
    'totNbPhasesPerNode_s','interRedGreenAfterPh_s');
load(fileQgParam,'minGreenTimeSec');
gt_init = initialGreenTimePerPhase_s;
minGT = minGreenTimeSec; %4sec

disp('rounding GTs if necessary');
gt_Opt = roundGTs10(fileData,fileNwk,fileQgParam,currPtId,loadGtsFromQg_bool,currGT_notFromQg,dbnIter,initPtId,dirName);


beforeRounding=0;
if loadGtsFromQg_bool
    print_gt_txt(fileOut_txt,fileData,fileNwk,fileQgParam,beforeRounding,gt_Opt);
else
    print_gt_txt_notFromQg(fileOut_txt,fileNwk,beforeRounding,gt_Opt);
end;
%{
expl:
gt_Opt = [23,4,4,12,16,4,21,...
        8,4,4,4,27,4,19,4,5,...
        31,28,13,...
        42,14,4,11,10,4,...
        4,17,13,22,4,...
        9,11,11,4,10,35,...
        63,22,...
        11,4,22,4,12,7,10,5,...
        4,35,4,4,16];
%}

if (sum(rem(gt_Opt,1))~=0) % decimal part of GTs
    error('GTs should now be integer values');
end;
% gt_Opt is of type double so test isinteger does not work

% check that all new greenTimes are greater than the minimal greenTime
if sum(gt_Opt < minGT)~=0
    error('there are green time smaller than the lower bound');
end;
% check that sum(gt) == available gt
for n=1:length(signalizedNode_s)
    currNodeID = signalizedNode_s(n);
    currPhases = find(phaseInfo(:,2)==currNodeID);
    if (abs(sum(gt_Opt(currPhases))-availableGreenTimes_s(n))>10^-12)
        error('sum gts inconsistent with available gts');
    end;
end;
% check that gt are integers
if (floor(gt_Opt)~= gt_Opt)
    error('should we set all green times to integer values?');
end;

% the phase Indexes (global indexes) that have new gt's (new control plan)
phDiff = find(gt_Opt~=gt_init);
% get their corresponding node indexes
newGT_nodeIDs = phaseInfo(phDiff,2)';
newGT_nodeIDs = unique(newGT_nodeIDs);

allGt_opt ={};
allGt_init ={};
allDiffStartT_opt={}; % the difference in the start time wrt to the initial plan
allStartT_opt={};
allStartT_init={};

% FOR A GIVEN NODE
if isempty(newGT_nodeIDs)
    disp('rounded gt == gt_init');
else
    for currNode_ID = newGT_nodeIDs
        
        % get the index of that node
        currNodeIndx = find(signalizedNode_s==currNode_ID);
        
        % get the set of phase IDs (global IDs) of this node
        curr_phase_globalIDs = find(phaseInfo(:,2)==currNode_ID);
        % get the local phase IDs that they correspond to.
        curr_phase_localIDs = phaseInfo(curr_phase_globalIDs,1);
        curr_fixedPhase_localIDs = fixedPhases_T_s{currNodeIndx}(1,:); % set oh phases that have a fixed gt
        curr_fixedPhase_gt = fixedPhases_T_s{currNodeIndx}(2,:); % their corresponding gt
        % total nb of phases for this node (both vble and fixed)
        currTotNbPhases = totNbPhasesPerNode_s(currNodeIndx);
        
        % since in aimsun no offsets are determined we set the offset of the node
        % (intitial start time of the whole control plan) to that of phase 1
        offset_phaseID = 1; %min(curr_phase_localIDs);
        
        % for all phases in the node (vble or fixed) update their initialStartTime and their greenTimeDuration
        % (for fixed phases we also shift their startTime)
        
        allGt_opt{currNodeIndx} = zeros(1,currTotNbPhases);
        allGt_opt{currNodeIndx}(curr_phase_localIDs) = gt_Opt(curr_phase_globalIDs);
        allGt_opt{currNodeIndx}(curr_fixedPhase_localIDs) = curr_fixedPhase_gt;
        
        allGt_init{currNodeIndx} = zeros(1,currTotNbPhases);
        allGt_init{currNodeIndx}(curr_phase_localIDs) = gt_init(curr_phase_globalIDs);
        allGt_init{currNodeIndx}(curr_fixedPhase_localIDs) = curr_fixedPhase_gt;
        
        % add the interRedGreen times
        allStartT_opt{currNodeIndx} = [0,cumsum(allGt_opt{currNodeIndx})] + cumsum([0,interRedGreenAfterPh_s{currNodeIndx}]);
        allStartT_init{currNodeIndx} = [0,cumsum(allGt_init{currNodeIndx})] + cumsum([0,interRedGreenAfterPh_s{currNodeIndx}]);
        
        allDiffStartT_opt{currNodeIndx} = allStartT_opt{currNodeIndx} - allStartT_init{currNodeIndx};
    end;
end;

% print the results to feed to aimsun - python script
printGT_pyth10(fileOut_py,newGT_nodeIDs,signalizedNode_s,totNbPhasesPerNode_s,allDiffStartT_opt,allGt_opt,dirName);
save(gtOptFile,'gt_Opt');

