function [x_out,fval,exitflag,output,lam_lagrange,grad_fmin,hessian_fmin]= fmincon_totTimeGamma_someMuExog_cf6_Phi_onlyQg_i(xk,loadX0_inF,qgParamsF,nwkF,outFilePrint,saveX0_outF,myOptionsOnlyQg)
% differs from fmincon_totTimeGamma_someMuExog_cf6_Phi_quad_i: only solves fmincon for qg model (no polynomial)

% xk contains the current green times
% loadX0: file with x_init

%%%%%
% Preliminaries
%%%%%
dbstop if error
dbstop if warning

% afaf are all of these vbles really necessary? exclude those that are only needed for calculating other ones
load(myOptionsOnlyQg,'myOptions');


load(qgParamsF,'minGreenTimeSec','nb','nbVbles','mu_index','last_index','nbPhases','nbSignNodes','nbSignQs');

load(nwkF,'gamma_s','serv_s','buffer_s','pij_s','mu_s','freeFlowCap','availableGreenTimes_s','fixedGreenTime','phase2Node',...
    'q2phases','cycleTimes_s','signalizedQs_s','signalizedQs2Nodes_s','condFact');

gamma = gamma_s;
serv = serv_s;
buffer = buffer_s;
pij = pij_s;
mu = mu_s;
cycleTime = cycleTimes_s;
signalizedQs = signalizedQs_s;
signalizedQs2Nodes = signalizedQs2Nodes_s;
availGreenTime = availableGreenTimes_s;
clear *_s;

[term_index,non_term_index,nbSucc,maxSucc,succAll,nbPred,maxPred,predAll] = prelim_storage_cf(serv, pij, buffer);

% knownX0 == 1
% xk is always a past iterate (ie derived as the solution of the optimization problem)
% so coherent qg parameters exist.
disp(['xInit for fmincon loaded from file: ',loadX0_inF]);
load(loadX0_inF,'x_out');
x_init = x_out;
clear x_out;

gt_indx = last_index:(last_index+nbPhases-1);
allMu = mu;

if (sum(xk ~= x_init(gt_indx))>0) error('incompatible green times'); end;


%%%%%
% linear  equality constraints.
%%%%%



%%%
% A eq
%%%

% normalized vsn:
% eqtns : sum(green times of phases) = availableGreenTime/cycleTime of n^th node
rows_1 = phase2Node;
cols_1 = last_index+(1:nbPhases)-1;
vals_1 = ones(1,nbPhases);
% normalized vsn:
% eqtns = 'mu - sum(green times of phases)*freeFlowCap == fixedGreenTime*freeFlowCap/cycleTime '
rows_2 = [];
cols_2 = [];
vals_2 = [];
disp('af vectorize this');
for j =1:length(signalizedQs)
    q = signalizedQs(j);
    
    currPhases = q2phases{q};
    currRow = nbSignNodes+j;
    
    rows_2 = [rows_2,currRow*ones(1,length(currPhases)+1)];
    
    muCol = mu_index+j-1;
    phaseCols = last_index+currPhases-1;
    cols_2 = [cols_2,muCol,phaseCols];
    
    currNodeIndx = phase2Node(currPhases);
    currNodeIndx = unique(currNodeIndx);
    if (length(currNodeIndx)~=1) error('all phases should be associated to the same node'); end;
    
    valT = freeFlowCap;
    vals_2 = [vals_2,1,-valT*ones(1,length(currPhases))];
end;

% Pif diag elements:
rows_3 = nbSignNodes + nbSignQs + (1:nb);
cols_3 = nb*4 + (1:nb);
vals_3 = ones(1,nb);
% Pif non terminal
rows_4 = [];
cols_4 = [];
vals_4 = [];
for currClass_succ = 1:maxSucc
    % only deal with non terminal queues
    currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ));
    % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
    
    % get the indexes of their successors
    currSucc_ind = succAll(currQs_1,1:currClass_succ);
    %disp('all of the above values should be calculated during the preprocessing phase');
    
    for currS_id = 1:currClass_succ
        % calc the linear indexes in order to access the element:
        % pij(currQs_1,currSucc_ind(:,currS_id))
        lin_indx = (currSucc_ind(:,currS_id)'-1)*nb+currQs_1;
        
        rows_4 = [rows_4,nbSignNodes + nbSignQs + currQs_1];
        cols_4 = [cols_4,nb*3+currSucc_ind(:,currS_id)'];
        vals_4 = [vals_4,-pij(lin_indx)];
    end;
    clear currQs_1 currSucc_ind lin_indx
end;


Aeq = sparse([rows_1,rows_2,rows_3,rows_4],...%,rows_5,rows_6,rows_7,rows_8,rows_9,rows_10],...
    [cols_1,cols_2,cols_3,cols_4],...%cols_5,cols_6,cols_7,cols_8,cols_9,cols_10],...
    [vals_1,vals_2,vals_3,vals_4],...%vals_5,vals_6,vals_7,vals_8,vals_9,vals_10],
    nbSignNodes+nbSignQs+nb*1,nbVbles);
Aeq = full(Aeq);
clear rows_* cols_* vals_*;

%Aeq(signalizedIndx,signalizedIndx) = incidenceMatrix*cycleTime/freeFlowCap;


%%%
% b eq
%%%

% normalized vsn:
% beq for eqtns : sum(green times of phases) = availableGreenTime/cycleTime of n^th node
rows_1 = 1:nbSignNodes;
vals_1 = availGreenTime./cycleTime;
% normalized vsn:
% beq for eqtns = 'mu - sum(green times of phases)*freeFlowCap == fixedGreenTime*freeFlowCap/cycleTime '
rows_2 = nbSignNodes + (1:nbSignQs);
vals_2 = fixedGreenTime(signalizedQs).*(freeFlowCap./cycleTime(signalizedQs2Nodes));
% beq for etqns  'mu == mu_init'
% Pif: beq:=0

beq = sparse([rows_1,rows_2],1,[vals_1,vals_2],nbSignNodes+nbSignQs+nb*1,1);
beq = full(beq);
clear rows_* cols_* vals_*;

% set lower bound
nonNullCte = 10^-10;
lb = sparse(1,[last_index+(1:nbPhases)-1,nb+(1:nb)],[minGreenTimeSec./cycleTime(phase2Node),nonNullCte*ones(1,1*nb)],1,nbVbles);
%lb = sparse(1,last_index+(1:nbPhases)-1,minGreenTimeSec./cycleTime(phase2Node),1,nbVbles);
% constr: phase <= gt_min/cycleT  : Its a constr on phase vbles not on the mus
%lb(first_index+nb*5+signalizedIndx-1) = (minGreenTimeSec*freeFlowCap)./cycleTime(signalizedIndx);
disp('we could eventually set quite a few of the bounds to -Inf. try out.');
% ie use a lower bound for pi and tg but not for the other vbles.
% lower bound for t_G (==mu*CycleTime/S, where S is the freeFlow capacity)
ub = [];

% there are inconsistencies between the green time durations of the aimsun GUI and those from getDuration() (used to create outJunction3.txt and thus get the initial mu calculation)
% we base ourselves on outJunction3.txt values.
null_beq = find(beq(nbSignNodes+(1:nbSignQs))==0);
indx1 = find(abs(Aeq(nbSignNodes+null_beq,:)*(x_init'))>myOptions{5})';
nonNull_beq = find(beq(nbSignNodes+(1:nbSignQs))~=0);
indx2 = find(abs(Aeq(nbSignNodes+nonNull_beq,:)*(x_init')-beq(nbSignNodes+nonNull_beq))>myOptions{5})';
diffLines = [];
if ~isempty(indx1)  diffLines = nbSignNodes+null_beq(indx1); end;
if ~isempty(indx2)  diffLines = [diffLines;nbSignNodes+nonNull_beq(indx2)]; end;

if ~isempty(diffLines)
    error('some fixedGreenTime values might need to be changed so that total green times are consistent with those of outJunction3.txt');
    % save('Aeq_adjustGreenTimes.mat'); Cf. differenceGreenTimes_adjust.m
end;

[mtmp,itmp]=max(abs(Aeq*(x_init')-beq))
clear indx1 indx2 diffLines nonNull_beq null_beq curr* tmp*;

%keyboard; %save('tmpX_init_fmincon.mat','x_init');


% cte jacobian values
nbTerm_t = length(term_index);
nbNonTerm_t = length(non_term_index);
%nbSignTerm_t = length(intersect(term_index,signalizedQs));
nbJacInit = 2*nb+2*nbTerm_t + nbNonTerm_t;
jac_xInit = zeros(1,nbJacInit); % memory allocation = zeros(1,nzmax)
jac_yInit = zeros(1,nbJacInit);
jac_valsInit = zeros(1,nbJacInit);
% diagonal elements lambda, Pif
jac_xInit(1:(2*nb)) = (nb*2) + (1:(2*nb));
jac_yInit(1:(2*nb)) = (nb*2) + (1:(2*nb));
jac_valsInit(1:(2*nb)) = 1;
% diagonal elements for term Qs: mut, muHat
jac_xInit(2*nb+(1:nbTerm_t)) = term_index;
jac_yInit(2*nb+(1:nbTerm_t)) = term_index;
jac_valsInit(2*nb+(1:nbTerm_t)) = 1;
jac_xInit(2*nb+nbTerm_t+(1:nbTerm_t)) = nb+term_index;
jac_yInit(2*nb+nbTerm_t+(1:nbTerm_t)) = nb+term_index;
jac_valsInit(2*nb+nbTerm_t+(1:nbTerm_t)) = 1;
% deriv muHat diag element
jac_xInit(2*nb+2*nbTerm_t+(1:nbNonTerm_t)) = nb+non_term_index;
jac_yInit(2*nb+2*nbTerm_t+(1:nbNonTerm_t)) = nb+non_term_index;
jac_valsInit(2*nb+2*nbTerm_t+(1:nbNonTerm_t)) = 1;

%{
% jacobian of mu_hat wrt mu for term and signalized stns (endogMu)
[currQs, currTermSignIndx, currSignInd] = intersect(term_index,signalizedQs);
jac_xInit(2*nb+2*nbTerm_t+nbNonTerm_t+(1:nbSignTerm_t)) = nb+currQs;
jac_yInit(2*nb+2*nbTerm_t+nbNonTerm_t+(1:nbSignTerm_t)) = nb*5+currSignInd;
jac_valsInit(2*nb+2*nbTerm_t+nbNonTerm_t+(1:nbSignTerm_t)) = -1;
clear currQs currTermSignIndx currSignInd nbSignTerm_t nbTerm_t
%}

% nb of non null elements in jacobian
sumSucc = sum(nbSucc);
sumPred = sum(nbPred);
%nbNonSourceQ = find(nbPred>0);
%nbSignTerm_t = length(intersect(term_index,signalizedQs));
nbNonNullJac = nb*3 + length(non_term_index)*(2+2) + sumSucc*2 + sumPred*1 + length(signalizedQs);
%nbNonNullJac = nb*3 + length(non_term_index)*(2+2) + sumSucc*2 + sumPred*1 + length(intersect(non_term_index,signalizedQs));
clear sumSucc sumPred;

 

%NOTE gradient is defined as the transposed of that of fsolve
% options=optimset('Diagnostics','on','Algorithm',myOptions{2},'GradObj','on','GradConstr','on','Display','iter','MaxFunEvals',100000000,'MaxIter',myOptions{3},...
%     'TolFun',myOptions{1},'TolCon',myOptions{5},'DerivativeCheck',myOptions{4}); %,'RelLineSrchBnd',myOptions{2}); %,'TolX',1e-36,'TolFun',1e-36);
options=optimset('Diagnostics','on','Algorithm',myOptions{2},'GradObj','on','GradConstr','on','Display','iter','MaxFunEvals',100000000,'MaxIter',myOptions{3},...
    'TolFun',myOptions{1},'TolCon',myOptions{5},'DerivativeCheck',myOptions{4}); %,'RelLineSrchBnd',myOptions{2}); %,'TolX',1e-36,'TolFun',1e-36);
% useful only for large-scale algo:
%'Hessian','on',
%'Algorithm','interior-point'
%'DerivativeCheck','on','Diagnostics','on',
%'Hessian','user-supplied','HessFcn',@objFctnHess,


x_init2 = x_init;
x_init2(nb+(1:nb)) = 1./x_init(nb+(1:nb));
tic
[x_out,fval,exitflag,output,lam_lagrange,grad_fmin,hessian_fmin] = fmincon(@objFctn,x_init2,[],[],Aeq,beq,lb,ub,@nonLinConstr,options);
fMinTime = toc;
rho = x_out(nb*2+(1:nb)).*x_out(nb+(1:nb));
x_out(nb+(1:nb)) = 1./x_out(nb+(1:nb));


if ( nargin>6 )
    print_fmincon_out(fval, exitflag, output, lam_lagrange, outFilePrint);
end;

save(saveX0_outF,'x_out','fMinTime');
xk_New = x_out(gt_indx);

if (exitflag <= 0) error('fmincon did not cv'); end;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
    function [z, z_grad] = objFctn(y0)
        %function [z, z_grad,z_hessian] = objFctn(y0)
        % returns a scalar.
        
        % TT := sum(Ni) / sum(gamma(1-pfull))
        
        rho = y0(nb*2+(1:nb)).*y0(nb+(1:nb));
        allNi_z2 = sum( rho.*(1./(1-rho)-( ((buffer+2).*(rho.^(buffer+1)))./(1-rho.^(buffer+2)) )));
        
        hasGammaInd_z2 = find(gamma>0);
        allG_z2 = sum(gamma(hasGammaInd_z2).*(1-y0(nb*3 + hasGammaInd_z2)));
                
        z = (allNi_z2/((condFact/60)*allG_z2));
        
        %gradient
        derivRho = 1./((1-rho).^2) - (((buffer+2).^2).*(rho.^(buffer+1)))./((1-(rho.^(buffer+2))).^2);
        jac_val_Q = (1/(condFact/60))*[(allNi_z2.*gamma(hasGammaInd_z2))/((allG_z2)^2),...
            (derivRho.*y0(nb*2+(1:nb)))/allG_z2,...
            (derivRho.*y0(nb+(1:nb)))/allG_z2];
        jac_indx_Q = [nb*3+hasGammaInd_z2,nb+(1:nb),nb*2+(1:nb)];
                
        % ensure that indexes of gradients are not overlapping
        z_grad = sparse(1,jac_indx_Q,jac_val_Q,1,nbVbles);
        
        % hessian
        % z_hessian=sparse(muHatIndx,muHatIndx,2./(y0(muHatIndx).^3),nbVbles,nbVbles);
        clear tmp* clear *z2 jac_val_* jac_indx_*;
        
    end %objFctn

% for interior-point algo the hessian must be calculated in a separate fctn: its not just the hessian of the objective fctn. Cf documentation.



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
    function [c, eq_out, grad_c, jacob] = nonLinConstr(s0)
        % OUTPUT:
        % c: compute nonlinear inequalities at x and their gradient:grad_c.
        % ceq : compute nonlinear equalities at x and their gradient: grad_ceq.
        
        %%%% coding of s0:
        % mu_t0 = s(nb+1,1:nb); unblocking rate of first blkd job.
        % mu_hat = s(nb+2,1:nb);
        % lamEFF = s(nb+3,1:nb); %INSTEAD OF lambdaPiQ
        % pj_full = s(nb+4,1:nb);
        % -- in linear equations:  pj_full_succ = s(nb+5,1:nb);
        % P(N>=1) = P(A=1)+P(B=1) -- overdetermined syst?
        %%%%
        nbEq = 4;
        eq_out=zeros(1,nbEq*nb);
        c = [];
        grad_c = [];
        %includes mut, muHat, lamPiQ, pjfull, eqtn:{P(N>=1)=P(A=1)+P(B=1)}
        % DOES NOT INCLUDE Pif, mu, phase variables
        
        % jacob initialization is done at the end (improves performance)
        jac_x = zeros(1,nbNonNullJac); % memory allocation = zeros(1,nzmax)
        jac_y = zeros(1,nbNonNullJac);
        jac_vals = zeros(1,nbNonNullJac);
        currJacInd = 0;
        
        % diagonal elements: calculated in preprocesing
        
        %%%%%
        % updating of syst vbles
        %%%%
        
        % usefull to evaluate eq_out without always distinguishing between muExog and muEndog.
        % thus we on ly need to distinguish when calculating jacobian
        allMu(signalizedQs) = s0(nb*5+(1:nbSignQs));
        rho = s0(nb*2+(1:nb)).*s0(nb+(1:nb));
        
        if ~isempty(term_index) % special cases of terminal stns
            % mut
            eq_out(term_index)=s0(term_index); % unblocking rate assumed null (CONVENTION, because theoretically it should be infinite)
            
            % mu_hat
            eq_out(nb+term_index)=s0(nb+term_index) - 1./allMu(term_index);
            
            % jacobians are calculated in the preprocessing
            
            [currQs, currTermSignIndx, currSignInd] = intersect(term_index,signalizedQs);
            currLeng = length(currQs);
            currAllVals = currJacInd+(1:currLeng);
            jac_x(currAllVals) = nb+currQs;
            jac_y(currAllVals) = nb*5+currSignInd;
            jac_vals(currAllVals) = 1./(s0(nb*5+currSignInd).^2);
            currJacInd = currJacInd+currLeng;
            clear currQs currTermSignIndx currSignInd;
            
        end;
        
        
        %disp(' Cf values ci-dessous, they should be calculated during the preprocessing phase');
        % mut
        % calculate mut equations vectorized depending on the nb of successors of each queue
        denomMut_class = zeros(1,nb);
        for currClass_succ = 1:maxSucc
            % only deal with non terminal queues
            currQs_1 = find(nbSucc==currClass_succ)';
            % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
            % currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ))';
            
            % get the indexes of their successors
            currSucc_ind = succAll(currQs_1,1:currClass_succ);
            %disp('all of the above values should be calculated during the preprocessing phase');
            denomMut_class1 = s0(nb*2+currSucc_ind).*s0(nb+currSucc_ind);
            
            
            % sum for a given queue across its successor queues
            % sum across ROWS
            if currClass_succ==1
                denomMut_class(currQs_1) = denomMut_class1;
            else
                denomMut_class(currQs_1) = sum(denomMut_class1,2);
            end;
            % full mut equation is calculated outside of this loop
            
            
            currLeng = length(currQs_1);
            for currS_id = 1:currClass_succ
                
                % jacobian : mut
                currAllVals = currJacInd+(1:currLeng);
                jac_x(currAllVals)=currQs_1;
                jac_y(currAllVals)=nb*2+currSucc_ind(:,currS_id);
                jac_vals(currAllVals) = -s0(nb+currSucc_ind(:,currS_id));
                currJacInd = currJacInd+currLeng;
                
                currAllVals = currJacInd+(1:currLeng);
                jac_x(currAllVals)=currQs_1;
                jac_y(currAllVals)=nb+currSucc_ind(:,currS_id);
                jac_vals(currAllVals) = -s0(nb*2+currSucc_ind(:,currS_id));
                currJacInd = currJacInd+currLeng;
            end;
            clear currQs_1 currSucc_ind denomMut_class1
        end;
        
        % mut eqtn
        eq_out(non_term_index)=s0(nb*2+non_term_index)./s0(non_term_index)-denomMut_class(non_term_index);
        
        currLeng = length(non_term_index);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = non_term_index;
        jac_y(currAllVals) = non_term_index;
        jac_vals(currAllVals) = -s0(nb*2+non_term_index)./(s0(non_term_index).^2);
        currJacInd = currJacInd+currLeng;
        
        currLeng = length(non_term_index);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = non_term_index;
        jac_y(currAllVals) = nb*2+non_term_index;
        jac_vals(currAllVals) = 1./s0(non_term_index);
        currJacInd = currJacInd+currLeng;
        
        % muHat
        % eqtn is vectorized for single server stations
        
        eq_out(nb+non_term_index)=s0(nb+non_term_index)-(1./allMu(non_term_index)+s0(nb*4+non_term_index)./s0(non_term_index));
        
        currLeng = length(non_term_index);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb+non_term_index;
        jac_y(currAllVals) = non_term_index;
        jac_vals(currAllVals) = s0(nb*4+non_term_index)./(s0(non_term_index).^2);
        currJacInd = currJacInd+currLeng;
        
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb+non_term_index;
        jac_y(currAllVals) = nb*4+non_term_index;
        jac_vals(currAllVals) = -1./s0(non_term_index);
        currJacInd = currJacInd+currLeng;
        
        % derivative of mu_hat wrt mu
        % non_term stns + endogMu
        [currQs, currNonTermSignIndx, currSignInd] = intersect(non_term_index,signalizedQs);
        currLeng = length(currQs);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb+currQs;
        jac_y(currAllVals) = nb*5+currSignInd;
        jac_vals(currAllVals) = 1./(s0(nb*5+currSignInd).^2);
        currJacInd = currJacInd+currLeng;
        
        % term stns + endogMu: % calculated in preprocessing
        
        % lambda: lamPiQ
        
        % source stns
        sourceQ = find(nbPred==0);
        eq_out(nb*2+sourceQ)=s0(nb*2+sourceQ) - gamma(sourceQ).*(1-s0(nb*3+sourceQ));
        
        %jacobian
        currLeng = length(sourceQ);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb*2+sourceQ;
        jac_y(currAllVals) = nb*3+sourceQ;
        jac_vals(currAllVals) = gamma(sourceQ);
        currJacInd = currJacInd+currLeng;
        
        % non-source stns
        nonSourceQ = find(nbPred>0);
        
        %disp(' Cf values (ci-dessous), they should be calculated during the preprocessing phase');
        for currClass_pred = 1:maxPred
            currQs_p = find(nbPred==currClass_pred)';
            % contains the q IDs of the non_term queues with nb of predecessors equal to currClass_pred
            % currQs_p = nonSourceQ(find(nbPred(nonSourceQ)==currClass_pred))';
            
            % get the indexes of their predecessors
            currPred_ind = predAll(currQs_p,1:currClass_pred);
            
            pij_pred = zeros(length(currQs_p),currClass_pred);
            for currP = 1:currClass_pred
                % get the linear index of the elements
                % pij(currPred_ind(:,currP),currQs_p)
                lin_indx = ((currQs_p - 1)*nb + currPred_ind(:,currP))';
                %disp('all of the above values should be calculated during the preprocessing phase');
                pij_pred(:,currP) = pij(lin_indx).*s0(nb*2+currPred_ind(:,currP));
                
                % part of the jacobian for lambda
                currLeng = length(currQs_p);
                currAllVals = currJacInd+(1:currLeng);
                jac_x(currAllVals) = (nb*2+currQs_p);
                jac_y(currAllVals) = nb*2+currPred_ind(:,currP);
                jac_vals(currAllVals) = -pij(lin_indx);
                currJacInd = currJacInd+currLeng;
            end;
            
            % lambda eqtn.
            currSum = (sum(pij_pred,2))';
            eq_out(nb*2+currQs_p) = s0(nb*2+currQs_p) - (gamma(currQs_p).*(1-s0(nb*3+currQs_p)) + currSum);
            
            currLeng = length(currQs_p);
            currAllVals = currJacInd+(1:currLeng);
            jac_x(currAllVals) = nb*2+currQs_p;
            jac_y(currAllVals) = nb*3+currQs_p;
            jac_vals(currAllVals) =  gamma(currQs_p);
            currJacInd = currJacInd+currLeng;
            
            clear currQs_p currPred_ind pij_pred currP lin_indx
        end;
        
        clear nonSourceQ
        
        
        % pjfull using closed form eqtns for loss models
        eq_out(nb*3+(1:nb))= s0(nb*3+(1:nb)) - ((1-rho).*(rho.^(buffer+1)))./(1-(rho.^(buffer+2)));
        
        % pjfull jacobian
        % diag element is in preprocessing
        currLeng = nb;
        
        derivRho = ((rho.^buffer).*(buffer+1-(buffer+2).*rho))./(1-(rho.^(buffer+2)))+((buffer+2).*(rho.^(2.*(buffer+1))).*(1-rho))./((1-(rho.^(buffer+2))).^2);
        %/(1-(rho.^(buffer+2)))).*(-1 + (buffer+2).*s0(nb*3+(1:nb))) + (s0(nb*3+(1:nb)).*(buffer+1))./rho;
        % wrt lambda
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb*3+(1:nb);
        jac_y(currAllVals) = nb*2+(1:nb);
        jac_vals(currAllVals) = -derivRho.*s0(nb+(1:nb));
        currJacInd = currJacInd+currLeng;
        % wrt muHat
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb*3+(1:nb);
        jac_y(currAllVals) = nb+(1:nb);
        jac_vals(currAllVals) = -(derivRho.*s0(nb*2+(1:nb)));
        currJacInd = currJacInd+currLeng;
        
        % eqtn {P(N>=1)=P(A=1)+P(B=1)}
        %disp('I believe P(N>=1) eqtn makes the system of eqtns overdetered, thus im not including it for now');
        %eq_out(nb*4+(1:nb))=  1- (1-rho)./(1-(rho.^(buffer+2))) - ( ( s0(nb*2+(1:nb))./(s0(nb*2+(1:nb))+s0(nb+(1:nb))) ).*( 1+ (allMu.*s0(nb*4+(1:nb)))./(s0(1:nb)+allMu.*s0(nb*4+(1:nb)).*(s0(nb*2+(1:nb))./(s0(nb*2+(1:nb))+s0(nb+(1:nb)))) ) );
        %afaf derivs
        
        % check that nbJacInit is the exact value :
        if currJacInd ~= length(jac_x) error('recalculate the nb of jacobian elements initialized or copy jac_x (code follows)'); end;
        
        % for fmincon we output the transpose of the jacobian
        % jacob = sparse(jac_y,jac_x,jac_vals,nbVbles,nbVbles-nb-nbPhases); % does not include mu not phase vbles
        jacob = sparse([jac_y,jac_yInit],[jac_x,jac_xInit],[jac_vals,jac_valsInit],nbVbles,nbEq*nb); % does not include mu, Pif, phase vbles
        
        % used for fmincon
        % UNLIKE fsolve:  fmincon defines jacob(i,j) is the partial derivative of equation(j) with respect to x(i)
        % thus need to transpose jacob_2dim;
        % jacob = jacob';
        
        %grad_ceq=[jacob_2dim;zeros(1,nbVbles)]; % adding derivatives wrt to z (they are all null)
        %max(abs(eq_out))
        %cond(jacob)
        %condest(jacob) % for sparse structure
        
        if length(jac_x) > nbNonNullJac error('recalculate the nb of jacobian elements'); end;
        clear jac_vals jac_x jac_y s0 currAllVals %j_tmp; %jacob_2dim;
    end % fctn nonLinConstr
end % fctn netwk1eqtn
