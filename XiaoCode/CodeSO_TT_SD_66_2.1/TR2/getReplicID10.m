function replicID = getReplicID10(dirName)

fileN = [dirName,'PyOutputs/currRepID.txt'];
fid=fopen(fileN,'r');
replicID_cell = textscan(fid,'%d');
fclose(fid);

if length(replicID_cell)~=1
    error('');
else
    replicID = replicID_cell{1};
end;