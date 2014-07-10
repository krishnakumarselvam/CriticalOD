function [TrialPoint] = OptimizeMetamodel(Betas,baseODMatrix,x_init,NUM_VEHICLES_TO_REMOVE,TopODIndices,AllowedReductionPercentage)
    
%Here I optimize the quadratic metamodel with the constraints that the
%demand for each OD pair is >=0
%Added a constraint that exactly the required number of vehicles are removed 
% I amassuming a quadratic metamodel of the form y = Constant + Beta_(2 to n+1)*X +
% Beta(n+1 to 2n+1)*X^2

%Here I add the queuing model to the metamodel

RHS = baseODMatrix(:,3);
PROBLEMDIMENSION = length(TopODIndices);

function [z] = obj(x0)
    z= Betas(1) + Betas(2 : PROBLEMDIMENSION+1)*x0 + Betas(PROBLEMDIMENSION + 2 : 2*PROBLEMDIMENSION+1)*x0.^2;
end

A   = ones(1,PROBLEMDIMENSION);
lb  = (1-AllowedReductionPercentage) * baseODMatrix(TopODIndices,3);
ub  =baseODMatrix(TopODIndices,3);
RHS = sum(ub)-NUM_VEHICLES_TO_REMOVE;
options=optimset('Diagnostics','off','Algorithm','interior-point','Display','off','MaxFunEvals',100000000,'MaxIter',100000000,...
    'TolFun',10e-7,'TolCon',10e-7);
CurrTrialPoint = fmincon(@obj,x_init',[],[],A,RHS,lb,ub,[],options);
TrialPoint=CurrTrialPoint';

end