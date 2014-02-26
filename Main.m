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
MAXITER = 50;

%Reading input data
baseODMatrix = textread('Inputs/ODpairs.txt');%Read from the Original OD pair data
AllowedReductionPercentage = 0.25;
NUM_VEHICLES_TO_REMOVE  = 40; % If this value is changed, the corresponding value in GenerateRandomTrialPoints also needs to be changed

Evaluated_Points = zeros(MAXITER,size(baseODMatrix,1));
Fsimvalues = zeros(MAXITER,1);

%%
CurrBeta  = zeros(1,size(baseODMatrix,1));

%%
for iter = 1:MAXITER

    %%

    %Find a trial point
    %-----------------------------------------
    [TrialPoint]=FindTrialPoint(iter,baseODMatrix,HOMEDIRECTORY,CurrBeta,Evaluated_Points,NUM_VEHICLES_TO_REMOVE);
    %-----------------------------------------
    Evaluated_Points(iter,:)=TrialPoint;

    %%

    %Evaluate trial point
    Fsimvalues(iter,1) = GetFsim(iter,HOMEDIRECTORY);
    
    %Update metamodel
    CurrBeta = UpdateMetamodel(Fsimvalues,Evaluated_Points);
    CurrBeta = CurrBeta';


    %%

    %Accept / Reject trial point


    %%

    %Mixture points


end

%%
%Post processing
plot(Fsimvalues,'-*')
