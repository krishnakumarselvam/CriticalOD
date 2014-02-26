function [CurrTrialPoint] = OptimizeMetamodel(Betas,baseODMatrix,x_init,NUM_VEHICLES_TO_REMOVE)
    
%Here I optimize the quadratic metamodel with the constraints that the
%demand for each OD pair is >=0

RHS = baseODMatrix(:,3);


function [z] = obj(x0)
    z=0;
    z = z+ Betas*x0;
end
A   = ones(1,size(baseODMatrix,1));
lb  = zeros(size(baseODMatrix,1),1);
ub  = baseODMatrix(:,3);
RHS = sum(ub)-NUM_VEHICLES_TO_REMOVE;
options=optimset('Diagnostics','off','Algorithm','interior-point','Display','off','MaxFunEvals',100000000,'MaxIter',100000000,...
    'TolFun',10e-7,'TolCon',10e-7);
CurrTrialPoint = fmincon(@obj,x_init',-A,-RHS,[],[],lb,ub,[],options);
CurrTrialPoint=CurrTrialPoint';
end