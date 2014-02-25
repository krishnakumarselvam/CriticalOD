function newAvgTTName = Create_Python_Input(iter,sz)
        
fid=fopen(['angParams.py'],'w');

filename = ['textFileName = \''C:\\Users\\Krishna\\Dropbox\\CriticalOD\\Aimsun\\' sz '\\Iter' num2str(iter) '.txt\'''];
fprintf(fid,['newAngName = \''C:\\Users\\Krishna\\Dropbox\\CriticalOD\\Aimsun\\' sz '.ang\''\n']);
fprintf(fid,['baseAngName = \''C:\\Users\\Krishna\\Dropbox\\CriticalOD\\Aimsun\\' sz '_base.ang\''\n']);
fprintf(fid,filename);

fclose(fid);
