function [TrialPoint] = OptimizeMetamodel(Betas,baseODMatrix,x_init,NUM_VEHICLES_TO_REMOVE,TopODIndices)
    
%Here I optimize the quadratic metamodel with the constraints that the
%demand for each OD pair is >=0
%Added a constraint that exactly the required number of vehicles are removed 

RHS = baseODMatrix(:,3);


function [z] = obj(x0)
    z=0;
    z = z+ Betas*x0;
end

PROBLEMDIMENSION = length(TopODIndices);

A   = ones(1,PROBLEMDIMENSION);
lb  = zeros(PROBLEMDIMENSION,1);
ub  = baseODMatrix(TopODIndices,3);
RHS = sum(ub)-NUM_VEHICLES_TO_REMOVE;
options=optimset('Diagnostics','off','Algorithm','interior-point','Display','off','MaxFunEvals',100000000,'MaxIter',100000000,...
    'TolFun',10e-7,'TolCon',10e-7);
CurrTrialPoint = fmincon(@obj,x_init',[],[],A,RHS,lb,ub,[],options);
TrialPoint=CurrTrialPoint';

%Testing module
% Testing = 1;
% if(Testing)
%     close all
%     plot(ub,'-*r')
%     hold on
%     plot(x_init,'-*g')
%     plot(CurrTrialPoint,'-*b')
% end
%     

end