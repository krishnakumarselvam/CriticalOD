function [Fsimvalue, Fsimvalue_Secondary] = GetFsim(iter,HOMEDIRECTORY)

%This function creates the aimsun file according to the changed OD Matrix
%as per the current trial point (details of which are in currTextFile). We
%then run the new aimsun file, and then return the average travel time
%resulting out of this.

%First create the necessary python input files, so that CreateModel.py can
%go about creating the aimsun file for this trial point

fid=fopen([HOMEDIRECTORY '\\PythonFiles\\angParams.py'],'w');

temp = strrep(HOMEDIRECTORY,'\','\\');

fprintf(fid,['baseAngName = \''' temp '\\Base\\AimsunBase.ang\''\n']);
fprintf(fid,['newAngName = \''' temp '\\AimsunFiles\\Aimsun_' num2str(iter) '.ang\''\n']);
fprintf(fid,['ResultFileName = \''' temp '\\Outputs\\Iter_' num2str(iter) '.txt\''\n']);
fprintf(fid,['ODFileName = \''' temp '\\TrialPoints\\Iter_' num2str(iter) '.txt\''']);
ResultFileName = [temp '\\Outputs\\Iter_' num2str(iter) '.txt'];

cd BatchFile
!CreateCommand
while ~exist(ResultFileName)
!RunCommand
end
cd (HOMEDIRECTORY)
filecontents = textread(ResultFileName);
Fsimvalue=filecontents(1,3);
Fsimvalue_Secondary = filecontents(1,1);
end