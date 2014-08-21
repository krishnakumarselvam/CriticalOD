function [xk_New, mk_New] = fmincon_totTimeGamma_someMuExog_cf6_Phi_quad30B_i(xk,loadX0_inF,qgParamsF,nwkF,delta,betaFile,outFilePrint,saveX0_outF,fctnTol,constrTol,phases2Remove)
% differs from fmincon_totTimeGamma_someMuExog_cf6_Phi_quad_i: poly only excludes one phase per intersection

% xk contains the current green times
% loadX0: file with the current green times and
% betaFile is the name of the mat file that contains the params of the model mk: a,b,C and alpha

disp('norm2 is used for trust-region constraint');

%%%%%
% Preliminaries
%%%%%
dbstop if error
dbstop if warning

% afaf are all of these vbles really necessary? exclude those that are only needed for calculating other ones
load(qgParamsF,'myOptions','minGreenTimeSec','nbPhases','nbSignNodes');

load(nwkF,'availableGreenTimes_s','cycleTimes_s','phase2Node');
cycleTime = cycleTimes_s;
availGreenTime = availableGreenTimes_s;
clear *_s;

% knownX0 == 1
% xk is always a past iterate (ie derived as the solution of the optimization problem)
% so coherent qg parameters exist.
disp(['xInit for fmincon loaded from file: ',loadX0_inF]);
load(loadX0_inF,'x_out');
x_init = x_out;
clear x_out;

nbVbles = length(x_init);
%if (sum(xk ~= x_init(gt_indx))>0) error('incompatible green times'); end;
%keyboard;

load(betaFile,'betaOpt','aIndx','bIndx','cIndx');
nbDigits = length(int2str(ceil(max(abs(betaOpt)))));
scaleF = 10^(nbDigits-1);
if scaleF > 1
    aBeta = betaOpt(aIndx)/scaleF;
    bBeta = betaOpt(bIndx)/scaleF;
    cBeta = betaOpt(cIndx)/scaleF;
    disp([' **** scaled down beta params, scaleF: ',int2str(scaleF)]);
else
    disp([' **** beta params not scaled, scaleF would have been: ',int2str(scaleF)]);
    aBeta = betaOpt(aIndx);
    bBeta = betaOpt(bIndx);
    cBeta = betaOpt(cIndx);
end;


%%%%%
% linear  equality constraints.
%%%%%



%%%
% A eq
%%%

% normalized vsn:
% eqtns : sum(green times of phases) = availableGreenTime/cycleTime of n^th node
rows_1 = phase2Node;
cols_1 = 1:nbPhases;
vals_1 = ones(1,nbPhases);


Aeq = sparse(rows_1,cols_1,vals_1,nbSignNodes,nbVbles);
Aeq = full(Aeq);
clear rows_* cols_* vals_*;

%Aeq(signalizedIndx,signalizedIndx) = incidenceMatrix*cycleTime/freeFlowCap;


%%%
% b eq
%%%

% normalized vsn:
% beq for eqtns : sum(green times of phases) = availableGreenTime/cycleTime of n^th node
beq = (availGreenTime./cycleTime)';
clear rows_* cols_* vals_*;

% set lower bound
lb = minGreenTimeSec./cycleTime(phase2Node);
%lb = sparse(1,last_index+(1:nbPhases)-1,minGreenTimeSec./cycleTime(phase2Node),1,nbVbles);
% constr: phase <= gt_min/cycleT  : Its a constr on phase vbles not on the mus
%lb(first_index+nb*5+signalizedIndx-1) = (minGreenTimeSec*freeFlowCap)./cycleTime(signalizedIndx);
disp('we could eventually set quite a few of the bounds to -Inf. try out.');
% ie use a lower bound for pi and tg but not for the other vbles.
% lower bound for t_G (==mu*CycleTime/S, where S is the freeFlow capacity)
ub = [];

% there are inconsistencies between the green time durations of the aimsun GUI and those from getDuration() (used to create outJunction3.txt and thus get the initial mu calculation)
% we base ourselves on outJunction3.txt values.

[mtmp,itmp]=max(abs(Aeq*(x_init')-beq))
clear indx1 indx2 diffLines nonNull_beq null_beq curr* tmp*;

%keyboard; %save('tmpX_init_fmincon.mat','x_init');


%NOTE gradient is defined as the transposed of that of fsolve
options=optimset('Diagnostics','on','Algorithm',myOptions{2},'GradObj','on','GradConstr','on','Display','iter','MaxFunEvals',100000000,'MaxIter',myOptions{3},...
    'TolFun',fctnTol,'TolCon',constrTol,'DerivativeCheck',myOptions{4}); %,'RelLineSrchBnd',myOptions{2}); %,'TolX',1e-36,'TolFun',1e-36);
%options=optimset('Diagnostics','on','Algorithm',myOptions{2},'GradObj','on','GradConstr','on','Display','iter','MaxFunEvals',100000000,'MaxIter',myOptions{3},...
%    'TolFun',myOptions{1},'TolCon',myOptions{5},'DerivativeCheck',myOptions{4}); %,'RelLineSrchBnd',myOptions{2}); %,'TolX',1e-36,'TolFun',1e-36);
% useful only for large-scale algo:
%'Hessian','on',
%'Algorithm','interior-point'
%'DerivativeCheck','on','Diagnostics','on',
%'Hessian','user-supplied','HessFcn',@objFctnHess,


tic
[x_out,fval,exitflag,output,lam_lagrange,grad_fmin,hessian_fmin] = fmincon(@objFctn,x_init,[],[],Aeq,beq,lb,ub,@nonLinConstr,options);
fMinTime = toc;


if ( nargin>6 )
    print_fmincon_out(fval, exitflag, output, lam_lagrange, outFilePrint);
end;

save(saveX0_outF,'x_out','fMinTime');
xk_New = x_out;
% mk_New might differ from fval if the betas were scaled
if scaleF > 1
    mk_New = fval*scaleF;
else
    mk_New = fval;
end;

if (exitflag <= 0) error('fmincon did not cv'); end;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
    function [z, z_grad] = objFctn(y0)
        %function [z, z_grad,z_hessian] = objFctn(y0)
        % returns a scalar.
        
        % xk represents only the decision variables {green splits?/?times of the phases}
        % get the indices of the used phases
        phaseIndx = 1:length(y0);
        phaseIndx(phases2Remove)=[];
        z = aBeta + y0(phaseIndx)*bBeta+(y0(phaseIndx).^2)*cBeta;
        z_grad = zeros(1,length(y0));
        z_grad(phaseIndx) = bBeta + 2.*(y0(phaseIndx)').*cBeta;
        
    end %objFctn

% for interior-point algo the hessian must be calculated in a separate fctn: its not just the hessian of the objective fctn. Cf documentation.



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
    function [c, eq_out, grad_c, jacob] = nonLinConstr(s0)
        % OUTPUT:
        % c: compute nonlinear inequalities at x and their gradient:grad_c.
        % ceq : compute nonlinear equalities at x and their gradient: grad_ceq.
        
        eq_out=[];
        jacob = [];

        %%%%%
        % TR constraint
        %%%%
        
        % bound on the step which is equal to: xTrial - xk
        sumT = sqrt(sum((s0-xk).^2));
        c = sumT - delta;
        grad_c = ((s0-xk)./sumT)';
        % are we sure that all other endogenous variables should be taken into account here?

    end % fctn nonLinConstr
end % fctn netwk1eqtn
