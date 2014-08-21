function printGT_pyth10(file_h,newGT_nodeIDs_h,signalizedNode_h,totNbPhasesPerNode_h,allDiffStartT_opt_h,allGt_opt_h,dirName)
%{
we do the same loop for file_h and for 'generic_new_gtContDico.py' (the generic one is created so that it can be imported by another python file, using a static file name (ie indep of current iteration)
 we print the values as dictionaries of dictionaries:
i.e. for each node we have a dictionary which we access as dico[nodeID] (lets call this dictionary nodeDico)
then for a given nodeDic we access a dictionary for every phase  as follows:
nodeDico = dico[nodeID];
phaseInfo_dico = nodeDico[phaseID];
and phaseInfo_dico[0] contains the difference in start time of that phase
    phaseInfo_dico[1] contains the new green time of that phase

python dico syntax:
aDico = {855: {2: [12,4], 3: [99,19]}, 1216: {7: [15,5]}}
%}

currIt = 1;
while currIt <= 2
    
    switch currIt
        case 1
            currFile = file_h;
        case 2
            currFile = [dirName,'PyInputs/generic_new_gtContDico.py'];
    end;
    
    
    if ~strcmp(currFile,'')
        fid=fopen(currFile,'w');
    else fid=1;
    end;
    
    fprintf(fid,'controlDico = {');
    
    if isempty(newGT_nodeIDs_h)
        % ie gt plan equals the intial gt plan; then controlDico vble is empty
        fprintf(fid,'}');
    else
        
        disp('we use internal phase IDs, which == the GUIs IDs-1');
        % FOR A GIVEN NODE
        for currNode_ID = newGT_nodeIDs_h
            fprintf(fid,'%i: {',currNode_ID);
            
            % get the index of that node
            currNodeIndx = find(signalizedNode_h==currNode_ID);
            
            currTotNbPhases = totNbPhasesPerNode_h(currNodeIndx);
            for k = 1:currTotNbPhases
                % some cases of phases that habe gt==0 are omitted
                % (aimsun GUI bug, comes probably from a phase that was created and then deleted but numbering of the phases is not consistent)
                if (allGt_opt_h{currNodeIndx}(k) > 0)
                    %disp('we use internal phase IDs, which == the GUIs IDs-1');
                    fprintf(fid,'%i: [',k-1); % k is the phase ID of the current phase
                    fprintf(fid,'%i,%i',allDiffStartT_opt_h{currNodeIndx}(k),allGt_opt_h{currNodeIndx}(k)); % new startTIme and greenTime of that phase
                    fprintf(fid,']');
                    if (k < currTotNbPhases) % its not the last phase
                        fprintf(fid,', ');
                    end;
                end;
            end;
            fprintf(fid,'}');
            if (currNode_ID ~= newGT_nodeIDs_h(end)) % its not the last node
                fprintf(fid,',\n');
            else
                fprintf(fid,'}');
            end;
        end;
    end;
    
    if ~strcmp(currFile,'')
        fclose(fid);
    end;
    
    currIt = currIt + 1;
end;