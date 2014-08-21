function addPoint(xStar,pointID,add2file,currReplicID)
% read most recent point (predictor) and combine it with the previous ones. 
% readDataF :   data file to read containing the used gt_Opt vector
% pointID :     id of the point to add (eg iteration)
% add2file :    file will all points obeserved, 1 pt per line, each line also has Dim+1 columns

fid=fopen(add2file,'a');
fprintf(fid,'%d ',pointID);
fprintf(fid,'%d ',currReplicID);
fprintf(fid,'%d ',xStar);
fprintf(fid,'\n');
fclose(fid);

ptDim = length(xStar);
currStr = '';
for j=1:ptDim
    currStr =[currStr,'%d '];
end;
% check that identifiers are unique
fid=fopen(add2file,'r');
dataObs = textscan(fid,['%d %d ',currStr]);%%ONLY READ FIRST TWO FROM FID. AFTER ',' APEND currStr 
fclose(fid);
[r,c] = size(unique([dataObs{1},dataObs{2}],'rows'));
if r~=length(dataObs{1})
    error('there are duplicate point IDs');
end;
