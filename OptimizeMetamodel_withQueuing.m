function [TrialPoint] = OptimizeMetamodel_withQueuing(Betas,baseODMatrix,x_init,NUM_VEHICLES_TO_REMOVE,TopODIndices,AllowedReductionPercentage)
    
%Here I optimize the quadratic metamodel with the constraints that the
%demand for each OD pair is >=0
%Added a constraint that exactly the required number of vehicles are removed 
% I amassuming a quadratic metamodel of the form y = Beta(1)*qg + Beta_(2 to n+1)*X +
% Beta(n+1 to 2n+1)*X^2

%Here I add the queuing model to the metamodel
keyboard
RHS = baseODMatrix(:,3);
DEMANDLENGTH = length(TopODIndices);
% %Objective function
% --------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------
function [z] = obj(x)
    z= Betas(1)*x(DEMANDLENGTH+1) + Betas(2 : DEMANDLENGTH+1)*x(1:DEMANDLENGTH) + Betas(DEMANDLENGTH + 2 : 2*DEMANDLENGTH+1)*x((1:DEMANDLENGTH)).^2;
end


% Adding non linear constraints
% --------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------
function [c,ceq] = nonlcon(x)
    
    % Setting lower and upper bounds
    lb  = (1-AllowedReductionPercentage) * baseODMatrix(TopODIndices,3);
    Ineq1= x(1:DEMANDLENGTH)-lb;
    ub  =baseODMatrix(TopODIndices,3);
    Ineq2= -x(1:DEMANDLENGTH)+ub;
    
    %Getting the equality constraint     
    Eq1 = sum(x(1:DEMANDLENGTH)) - (sum(ub)-NUM_VEHICLES_TO_REMOVE);
    
    %Adding the queuing model
    Eq2 = sum(x(1:DEMANDLENGTH))-x(DEMANDLENGTH+1)^2;

    c = [Ineq1;Ineq2];
    ceq = [Eq1; Eq2];
end

% Getting an initial point
% --------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------
InitialPoint = [x_init sum(x_init)]';


options=optimset('Diagnostics','off','Algorithm','interior-point','Display','off','MaxFunEvals',100000000,'MaxIter',100000000,...
    'TolFun',10e-7,'TolCon',10e-7);
CurrTrialPoint = fmincon(@obj,InitialPoint ,[],[],[],[],[],[],@nonlcon,options);
TrialPoint=CurrTrialPoint(1:DEMANDLENGTH)';

end