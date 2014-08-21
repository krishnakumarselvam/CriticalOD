function newAvgTTName = angParams2python10(currPtId,dbnIter,initPtId,dirName)
        
fid=fopen([dirName,'PyInputs/angParams.py'],'w');
if fid<=0 error(''); end;

dirName2 = '';
dirName3 = '';
for k=1:length(dirName)
    if strcmp(dirName(k),'/')
        dirName2 = [dirName2,'\'];
        dirName3 = [dirName3,'\'];
    else
        dirName2 = [dirName2,dirName(k)];
        dirName3 = [dirName3,dirName(k)];
    end;
end;

fprintf(fid,'newAngName = \''%sAimsun\\CentreV5_getAvgTT_NoReplic_%i_%i_%i.ang\''\n',dirName2,initPtId,dbnIter,currPtId);
fprintf(fid,'newAvgTTFile = \''%sPyOutputs\\TTPerReplic_%i_%i_%i.txt\''\n',dirName2,initPtId,dbnIter,currPtId);
fclose(fid);

newAvgTTName = [dirName3,'PyOutputs\TTPerReplic_',int2str(initPtId),'_',int2str(dbnIter),'_',int2str(currPtId),'.txt'];
% used by matlab to know which file to read the avgTT's from.
