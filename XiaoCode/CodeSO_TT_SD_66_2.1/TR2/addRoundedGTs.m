function addRoundedGTs(loadFile,add2file,currReplicID,currPtId) %(,fileData,fileNwk,fileQgParam,)
% appends into the file the iteration ID, replication ID and then the rounded gt plan

if isnumeric(loadFile) % the gtRd value has been passed directly and not through a file name
    gt_Opt = loadFile;
else
    load(loadFile,'gt_Opt');
end;

ptDim = length(gt_Opt);
nbCol = ptDim+2; % iterID replicID currGT

% check file format
if exist(add2file) %ie its not the very first iteration
    currStr = [];
    for j = 1:nbCol
        currStr = [currStr,' %d'];
    end;
    
    fid=fopen(add2file,'r');
    dataObs = textscan(fid,currStr);
    fclose(fid);
    
    dataObs = cell2mat(dataObs);
    [nbObs,c] = size(dataObs);
    if c~=nbCol error('inconsistent size'); end;
end;


% 2. print into file
fid=fopen(add2file,'a');
fprintf(fid,'%d ',currPtId);
fprintf(fid,'%d ',currReplicID);
fprintf(fid,'%d ',gt_Opt);
fprintf(fid,'\n ');
fclose(fid);
