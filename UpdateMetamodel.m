function [ betaOpt ] = UpdateMetamodel(Fsimvalues,Evaluated_Points)
%This function fits a quadratic metamodel with Fsimvalues and Evaluated
%Points

%!! Here we don't weight the regression errors based on the trial point, and
%we need to correct for it at some point !!

nbVblesBeta = size(Evaluated_Points,2);
beta0 = zeros(nbVblesBeta,1);
weightAugData = 10e-2;
C = [Evaluated_Points;weightAugData*eye(nbVblesBeta)];
d = [Fsimvalues; weightAugData*zeros(nbVblesBeta,1)];
options=optimset('Diagnostics','on','LargeScale','off','MaxFunEvals',100000000,'MaxIter',100000000,...
        'TolFun',1e-7,'DerivativeCheck','off'); 
[betaOpt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,[],[],[],[],[],[],beta0,options);


end

