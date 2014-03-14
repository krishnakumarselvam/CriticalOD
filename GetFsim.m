function [Fsimvalue] = GetFsim(iter,HOMEDIRECTORY)

%This function creates the aimsun file according to the changed OD Matrix
%as per the current trial point (details of which are in currTextFile). We
%then run the new aimsun file, and then return the average travel time
%resulting out of this.

%First create the necessary python input files, so that CreateModel.py can
%go about creating the aimsun file for this trial point

fid=fopen([HOMEDIRECTORY '\\PythonFiles\\angParams.py'],'w');

fprintf(fid,['baseAngName = \''C:\\Users\\Krishna\\Dropbox\\CriticalOD\\Base\\AimsunBase.ang\''\n']);
fprintf(fid,['newAngName = \''C:\\Users\\Krishna\\Dropbox\\CriticalOD\\AimsunFiles\\Aimsun_' num2str(iter) '.ang\''\n']);
fprintf(fid,['ResultFileName = \''C:\\Users\\Krishna\\Dropbox\\CriticalOD\\Outputs\\Iter_' num2str(iter) '.txt\''\n']);
fprintf(fid,['ODFileName = \''C:\\Users\\Krishna\\Dropbox\\CriticalOD\\TrialPoints\\Iter_' num2str(iter) '.txt\''']);

cd 'C:/Program Files/TSS-Transport Simulation Systems/AIMSUN 6.1/'
!aconsole.exe -script C:/Users/Krishna/Dropbox/CriticalOD/PythonFiles/Create_Model.py 

ResultFileName = ['C:\\Users\\Krishna\\Dropbox\\CriticalOD\\Outputs\\Iter_' num2str(iter) '.txt'];
while ~exist(ResultFileName)
!aconsole.exe -script C:/Users/Krishna/Dropbox/CriticalOD/PythonFiles/simRepAll_newMdl_newReplic.py 
end
cd (HOMEDIRECTORY)
filecontents = textread(ResultFileName);
Fsimvalue=filecontents(1,3);
end