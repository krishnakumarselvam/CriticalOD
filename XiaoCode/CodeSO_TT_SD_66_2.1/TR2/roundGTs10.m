function newGT = roundGTs10(fileData,fileNwk,fileQg,currPtId,loadGtsFromQg_bool,currGT_notFromQg,dbnIter,initPtId,dirName)
% currGT_notFromQg is the equivalent of x_out(last_index:(last_index+nbPhases-1)). It will be multiplied by cycleTimes_s(phase2Node) to be converted to actual green times

load(fileNwk,'cycleTimes_s','phase2Node','availableGreenTimes_s');
nbPhases = length(phase2Node);

if loadGtsFromQg_bool
    load(fileData,'x_out');
    load(fileQg,'last_index');
    initGT = x_out(last_index:(last_index+nbPhases-1)).*cycleTimes_s(phase2Node);
else
    initGT = currGT_notFromQg.*cycleTimes_s(phase2Node);
end;
    

newGT = round(initGT);

% round all GTs to the nearest integer
% then check wether the sums of GTs for a given node == availableGTs
% if not round accordingly (such that rounding per phase is minimized)
for currNode = 1:length(unique(phase2Node))
    % for a given node;
    currPhases = find(phase2Node==currNode);
    initGT_p = initGT(currPhases);
    newGT_p = newGT(currPhases);
    
    while (sum(newGT_p)~= availableGreenTimes_s(currNode))
        diffSum = sum(newGT_p) - availableGreenTimes_s(currNode);
        diffGTs = newGT_p-initGT_p; 
        if diffSum>0
            % need to reduce the green times
            % find the GT that has been rounded upwards but is closest to 0.5, so that rounding it downards is minimal
            candidates = find(diffGTs>0);
            % make sure we respect lower bounds
            ommitGT = find(newGT_p(candidates)==4);
            if ~isempty(ommitGT)
                %indxOmmit = find(candidates==ommitGT);
                candidates(ommitGT) = [];
            end;
            [mTmp, indxCand] = min(initGT_p(candidates)/10); %first decimal value
            GTchosen_indx = candidates(indxCand);
            newGT_p(GTchosen_indx) = floor(initGT_p(GTchosen_indx));
        else
            % need to increase the green times
            % find the GT that has been rounded downwards but is closest to 0.5, so that rounding it upwards is minimal
            candidates = find(diffGTs<0);
            % rounding upwards cannot violate lower bounds            
            [mTmp, indxCand] = max(initGT_p(candidates)/10); %first decimal value
            GTchosen_indx = candidates(indxCand);
            newGT_p(GTchosen_indx) = ceil(initGT_p(GTchosen_indx));
        end;        
    end;
    newGT(currPhases) = newGT_p;
end;

h=figure;
stem(initGT);
hold on;
stem(newGT,'r');
title('gtOpt vs rounded gts (red)');
figName = [dirName,'Figs/gtOpt_roundedVsOrig_',int2str(initPtId),'_',int2str(dbnIter),'_',int2str(currPtId)];
saveas(h,figName);
close(h);