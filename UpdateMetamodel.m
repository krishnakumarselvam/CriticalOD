function [ betaOpt ] = UpdateMetamodel(Fsimvalues,Evaluated_Points,LastAcceptedPoint,iter)

%This function fits a quadratic metamodel with Fsimvalues and Evaluated
%Points

%Removing the zero values
Evaluated_Points = Evaluated_Points(1:iter,:);
Fsimvalues = Fsimvalues(1:iter,:);

%Getting the weights for the different points
nbVblesBeta = size(Evaluated_Points,2);
beta0 = zeros(nbVblesBeta,1);
weightAugData = 10e-2;
xkCenter = Evaluated_Points(LastAcceptedPoint,:);
distAll = Evaluated_Points-repmat(xkCenter,iter,1);
norm_tmp = sqrt(sum((distAll').^2));
weightAll = 1./(1+norm_tmp)';

%Weighting the x and the y values
for i=1:nbVblesBeta
    Evaluated_Points(:,i) = weightAll.*Evaluated_Points(:,i);
end
Fsimvalues = Fsimvalues.*weightAll;

%Need to use the weights to fit the metamodel

C = [Evaluated_Points;weightAugData*eye(nbVblesBeta)];
d = [Fsimvalues; weightAugData*zeros(nbVblesBeta,1)];
options=optimset('Diagnostics','off','LargeScale','off','MaxFunEvals',100000000,'MaxIter',100000000,...
        'TolFun',1e-7,'DerivativeCheck','off'); 
[betaOpt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,[],[],[],[],[],[],beta0,options);


end

