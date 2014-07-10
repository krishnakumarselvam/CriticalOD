function [ betaOpt ] = UpdateMetamodel(Fsimvalues,Evaluated_Points,LastAcceptedPoint,iter)

%This function fits a quadratic metamodel with Fsimvalues and Evaluated
%Points

%Removing the zero values
Evaluated_Points = Evaluated_Points(1:iter,:);
Fsimvalues = Fsimvalues(1:iter,:);

%Getting the weights for the different points
PROBLEMDIMENSION = size(Evaluated_Points,2);
NUM_EVAL_POINTS = size(Evaluated_Points,1)

%Adding the constant and the squared terms
NEW_DiMENSION = 1 + 2*PROBLEMDIMENSION;
Expanded_Evaluated_Points = zeros(NUM_EVAL_POINTS,NEW_DiMENSION);
New_Ones = ones(NUM_EVAL_POINTS,1);
Expanded_Evaluated_Points = [New_Ones, Evaluated_Points, Evaluated_Points.^2 ]

beta0 = zeros(NEW_DiMENSION,1);
weightAugData = 10e-2;
xkCenter = Expanded_Evaluated_Points(LastAcceptedPoint,:);
distAll = Expanded_Evaluated_Points-repmat(xkCenter,iter,1);
norm_tmp = sqrt(sum((distAll').^2));
weightAll = 1./(1+norm_tmp)';


%Weighting the x and the y values
for i=1:NEW_DiMENSION
    Expanded_Evaluated_Points(:,i) = weightAll.*Expanded_Evaluated_Points(:,i);
end
Fsimvalues = Fsimvalues.*weightAll;

%Need to use the weights to fit the metamodel

C = [Expanded_Evaluated_Points;weightAugData*eye(NEW_DiMENSION)];
d = [Fsimvalues; weightAugData*zeros(NEW_DiMENSION,1)];
options=optimset('Diagnostics','off','LargeScale','off','MaxFunEvals',100000000,'MaxIter',100000000,...
        'TolFun',1e-7,'DerivativeCheck','off'); 
[betaOpt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,[],[],[],[],[],[],beta0,options);


end

