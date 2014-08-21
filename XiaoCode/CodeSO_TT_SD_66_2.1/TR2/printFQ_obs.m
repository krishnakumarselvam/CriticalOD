function printFQ_obs(add2file,currPtId,currFQ,fQ_ETT,fQ_SD)
% print the observed values to the output file 

fid=fopen(add2file,'a');
fprintf(fid,'%d %5.2f %5.2f %5.2f\n',currPtId,currFQ,fQ_ETT,fQ_SD);
fclose(fid);

% check that identifiers are unique
fid=fopen(add2file,'r');
dataObs = textscan(fid,'%d %f %f %f');
fclose(fid);
if (length(unique(dataObs{1}))~=length(dataObs{1}))
    error('there are duplicate point IDs');
end;