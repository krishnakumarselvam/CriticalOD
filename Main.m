%This is the main document for the Critical OD pair problem

%Description
%The objective here is the following: Given that we can remove x vehicles
%from different ODs in the network, which OD pair is the best candidate from which we could
%remove these x vehicles to make maximum impact on congestion

%%
clc
clear
delete('Outputs/*')
delete('AimsunFiles/*')
delete('TrialPoints/*')


%%

%Initialize parameters of SO algorithm
HOMEDIRECTORY = pwd;
MAXITER = 100;

%Reading input data
baseODMatrix = textread('Inputs/ODpairs.txt');%Read from the Original OD pair data
AllowedReductionPercentage = 0.25;
numNewStations=20;

Evaluated_Points = zeros(MAXITER,numNewStations);
Fsimvalues = zeros(MAXITER,1);

%%
for iter = 1:MAXITER

    %%

    %Find a trial point
    %-----------------------------------------
    [ChangedODMatrix,rowIDs,currTextFilename]=FindTrialPoint(iter,baseODMatrix,AllowedReductionPercentage,numNewStations,HOMEDIRECTORY);
    %-----------------------------------------
    Evaluated_Points(iter,:)=rowIDs';

    %%

    %Evaluate trial point and update metamodel
    Fsimvalues(iter,1) = GetFsim(iter,HOMEDIRECTORY);


    %%

    %Accept / Reject trial point


    %%

    %Mixture points


end

%%
%Post processing
plot(Fsimvalues,'-*')
