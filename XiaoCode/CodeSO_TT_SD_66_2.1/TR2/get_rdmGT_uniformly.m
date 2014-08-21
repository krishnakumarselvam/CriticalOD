function gt_plan = get_rdmGT_uniformly(fileNwk)
% generate a rdm gt plan

rand('state', sum(100*clock)); % Initialize rand to a different state each time:

load(fileNwk,'cycleTimes_s','phase2Node','availableGreenTimes_s','minGreenTimeSec');
nbPhases = length(phase2Node);
nbNodes = length(unique(phase2Node));

gt_plan = zeros(1,nbPhases);

for currNode = 1:nbNodes
    % for a given node
    
    currPhases = find(phase2Node==currNode);
    b_n = availableGreenTimes_s(currNode)/cycleTimes_s(currNode);
    nbPhases_n = length(find(phase2Node==currNode));
    g_min_n = minGreenTimeSec./cycleTimes_s(currNode);
    
    [x,v] = randfixedsum(nbPhases_n ,1,b_n,g_min_n,b_n);
        
    gt_plan(currPhases) = x;
end;

% verify that 1) sum(gt)==b_n 2) gt>=g_min
for currNode = 1:nbNodes
    currPhases = find(phase2Node==currNode);    
    b_n = availableGreenTimes_s(currNode)/cycleTimes_s(currNode);
    g_min_n = minGreenTimeSec./cycleTimes_s(currNode);

    if ~isempty(find(gt_plan(currPhases)<g_min_n)) error('gts do not satisfy lower bound constraints'); end;
    if (abs(sum(gt_plan(currPhases))-b_n)>10^-7) error('incorrect sum gts'); end;
end;
