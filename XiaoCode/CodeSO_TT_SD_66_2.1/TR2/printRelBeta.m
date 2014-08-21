function printRelBeta(currRelBeta,add2file,currPtId)
% print the observed values to the output file 

fid=fopen(add2file,'a');
fprintf(fid,'%d %2.5f\n',currPtId,currRelBeta);
fclose(fid);

% check that identifiers are unique
fid=fopen(add2file,'r');
dataObs = textscan(fid,'%d %f');
fclose(fid);
if (length(unique(dataObs{1}))~=length(dataObs{1}))
    error('there are duplicate point IDs');
end;