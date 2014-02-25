clc
clear
sz='mixed';
% delete('Large/*.txt')
% delete('Small/*.txt')
% delete('Mixed/*.txt')
for i=1:1
    Create_Python_Input(i,sz);
    cd 'C:/Program Files/TSS-Transport Simulation Systems/AIMSUN 6.1/'
    !aconsole.exe -script C:/Users/Krishna/Dropbox/CriticalOD/Aimsun/Create_Model.py 
%    pause(10)
%    !aconsole.exe -script C:/Users/Krishna/Dropbox/CriticalOD/Aimsun/simRepAll_newMdl_newReplic.py 
    cd 'C:/Users/Krishna/Dropbox/NEW_SO/Results/Comparison';
end
a
%Send_Notification