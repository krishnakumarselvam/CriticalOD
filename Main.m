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

%Initialize constants of SO algorithm
HOMEDIRECTORY = pwd;
MAXITER = 30;

%Initializing Parameters
ETA = .001;
BETATOL = 10^-1; % threshold for relative change in the norm of metamodel parameters (beta)

%Reading input data
baseODMatrix = textread('Inputs/ODpairs.txt');%Read from the Original OD pair data
AllowedReductionPercentage = 0.25;
NUM_VEHICLES_TO_REMOVE  = 40; % If this value is changed, the corresponding value in GenerateRandomTrialPoints also needs to be changed
load ('Inputs/RandomTrialPoints', 'TopODIndices');
PROBLEMDIMENSION = length(TopODIndices);

Evaluated_Points = zeros(MAXITER,PROBLEMDIMENSION);
Fsimvalues = zeros(MAXITER,1);
AcceptedFsimValues = zeros (MAXITER,1);
AcceptedInSO =  zeros(MAXITER,1);
LastAcceptedPoint=1;
IsMixturePoint = zeros(MAXITER,1);


%%
CurrBeta  = zeros(1,PROBLEMDIMENSION);

%%
for iter = 1:MAXITER

    %%

    %Find a trial point
    %-----------------------------------------
    [TrialPoint]=FindTrialPoint(iter,baseODMatrix,HOMEDIRECTORY,CurrBeta,Evaluated_Points,NUM_VEHICLES_TO_REMOVE,TopODIndices,AllowedReductionPercentage);
    %-----------------------------------------
    Evaluated_Points(iter,:)=TrialPoint;
    

    %%

    %Evaluate trial point
    Fsimvalues(iter,1) = GetFsim(iter,HOMEDIRECTORY);
  

    %%

    %Accept / Reject trial point
    IsCurrentPointAccepted = AcceptanceStep(iter,Fsimvalues,ETA,Evaluated_Points,LastAcceptedPoint,CurrBeta');
    AcceptedInSO(iter,1)= IsCurrentPointAccepted;
    if(IsCurrentPointAccepted)
        LastAcceptedPoint = iter;
    end
    %%
    %Update metamodel
    OldBeta = CurrBeta;
    CurrBeta = UpdateMetamodel(Fsimvalues,Evaluated_Points,LastAcceptedPoint,iter);
    CurrBeta = CurrBeta';


    %%

    %Mixture points
    DoWeNeedMixturePoint = EvaluateChangeinBeta(OldBeta,CurrBeta,BETATOL);
    if(DoWeNeedMixturePoint && iter+1 <= MAXITER)
        AcceptedFsimValues(iter) = Fsimvalues(LastAcceptedPoint);
        iter = iter+1;
        [TrialPoint] = FindMixturePoint(iter,baseODMatrix,HOMEDIRECTORY,TopODIndices);
        Evaluated_Points(iter,:)=TrialPoint;
        Fsimvalues(iter,1) = GetFsim(iter,HOMEDIRECTORY);
        OldBeta = CurrBeta;
        CurrBeta = UpdateMetamodel(Fsimvalues,Evaluated_Points,LastAcceptedPoint,iter);
        CurrBeta = CurrBeta';
        IsMixturePoint(iter) = 1;
    end
    
  AcceptedFsimValues(iter) = Fsimvalues(LastAcceptedPoint);
  PlotAlgorithmProgression(Fsimvalues,AcceptedInSO,iter,AcceptedFsimValues,IsMixturePoint);
 

end
%print -djpeg100 MatlabPlot.jpg


