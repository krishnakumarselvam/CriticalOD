function currAvg = addObservtns3(readF,pointID,add2file,runAimsun,currReplicID)
% read most recent observation and combine it with the previous ones.
% readF :       file to read containing 4 columns: 3rd column contains the average travel time
% pointID :     id of the point to add (eg iteration)
% add2file :    file that has 2 columns: the first is the id of the point (eg iteration); the second is the average travel time at that point

if isnumeric(readF) % the fSim value has been passed directly and not through a file name
    currAvg = readF;
else
    
    if runAimsun ~= 0
        currFile = readF;
    else
        currFile = 'C:\MIT\Lab\KanchanaNanduri\CodeSO_TT_FC\Runs3\PyOutputs\TTPerReplic_dummy.txt';
    end;
    
    fid=fopen(currFile,'r');
    % columns: total travel time; tot nb of vehicles that exited the nwk (after the warmup-time); average travel time; replication ID
    dataIn = textscan(fid,'%f %f %f %d');
    fclose(fid);
    
    currTTs = dataIn{3};
    if runAimsun == 0
        % generate random observations
        currTTs = .1*rand(6,1).*currTTs;
    end;
    
    currAvg = mean(currTTs);
end;

fid=fopen(add2file,'a');
fprintf(fid,'%d %d %5.5f\n',pointID,currReplicID,currAvg);
fclose(fid);

% check that identifiers are unique
fid=fopen(add2file,'r');
dataObs = textscan(fid,'%d %d %f');
fclose(fid);
[r,c] = size(unique([dataObs{1},dataObs{2}],'rows'));
if r~=length(dataObs{1})
    error('there are duplicates of (ptID,replic)');
end;
