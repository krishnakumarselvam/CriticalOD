function [x_out,fval,exitflag,output,jacobian,extra_vbles_out] = nwk_fsolve_muHat_cf5(gamma, lambda_flowCvtn, mu, serv, buffer, pij,outFile,scenarios,myOptions,knownX0,loadX0,iter_rdm)
% differs from nwk_fsolve_muHat_cf4.m : rho = lamEff/muHat instead of lamEff/(muHat*(1-pifull))
% INPUTS
% gamma             vector of external arrival rates (of length number of stations)
% lambda_flowCvtn   vector of theoretical flow cvtn arrival rates (need to initialize lambda0)
% mu                vector of service rates
% serv              vector of number of servers per station
% buffers           vector of buffer sizes per stn
% pij               matrix of routing probabilities
% outFile           tabulates  fsolve algo exitflag and fvalues.

if ~isempty(find(serv~=1)) error('this code is intended for single server networks only'); end;

%%%%%
% Preliminaries
%%%%%
dbstop if error
dbstop if warning


nb = length(buffer); %nb of nodes

[term_index,non_term_index,nbSucc,maxSucc,succAll,nbPred,maxPred,predAll] = prelim_storage_cf(serv, pij, buffer);

nb_extra_unknowns = 5; % does not include mu nor phase vbles.
nbVbles = nb*nb_extra_unknowns;

x_init=zeros(1,nbVbles);

if (knownX0 == 0)
    disp('calculating x0');
    % in this initialization mu is considered exogenous. i.e we initialize vbles to satisfy state eqtns but decision vbles play no role in the initialization.
    x_init = x0Init_storage_d_lin_newForm_muHat_cf(scenarios, serv,pij,mu,lambda_flowCvtn,nb,nb_extra_unknowns,term_index,non_term_index,maxSucc,nbSucc,succAll,buffer);
    disp('x0 calculation done');
    %keyboard; %save('x0_LIN_simLo.mat','x0_pi','x0_param');
    if strcmp(scenarios{1},'piRdm')
        %if  greenTimeScenario ~= 0
        %    filePi = ['Jacobian/x0_LIN_PiRdm_greenTime',int2str(greenTimeScenario),'_iter',int2str(iterRdm),'.mat'];
        %else
        filePi = ['Jacobian/x0_PiRdm_iter',int2str(iter_rdm),'.mat'];
        %end;
        save(filePi,'x_init');
    end;
else
    disp('make sure that we have loaded muHat and not 1/muHat');
    load(loadX0);
    %%%%%
    %% 1/muHat
    %x_init(nb + (1:nb)) = 1./x_init(nb + (1:nb));
    if length(x_init)~=nbVbles
        error('x0 is now loaded as a vector no longer as x0_param, and x0_pi so that we avoid recopying it');
    end;
    keyboard;
end;


% nb of non null elements in jacobian
sumSucc = sum(nbSucc);
sumPred = sum(nbPred);
%nbNonSourceQ = find(nbPred>0);
nbNonNullJac = nb*2 + length(non_term_index)*(3+2) + sumSucc*2 + sumPred*1 + length(find(nbPred>0))*1;
clear sumSucc sumPred;


%%%%%%
%%%%%%
% calculate all the jacobian values that are constant.
%%%%%%
%%%%%%

nbTerm_t = length(term_index);
%sourceQ = find(nbPred==0);
%nbSourceQ = length(sourceQ);
nbJacInit = 3*nb+2*nbTerm_t+sum(nbSucc);
jac_xInit = zeros(1,nbJacInit); % memory allocation = zeros(1,nzmax)
jac_yInit = zeros(1,nbJacInit);
jac_valsInit = zeros(1,nbJacInit);

% diagonal elements lambda, pjfull, Pif
jac_xInit(1:(3*nb)) = (nb*2) + (1:(3*nb));
jac_yInit(1:(3*nb)) = (nb*2) + (1:(3*nb));
jac_valsInit(1:(3*nb)) = 1;
% diagonal elements for term Qs: mut, muHat
jac_xInit(3*nb+(1:nbTerm_t)) = term_index;
jac_yInit(3*nb+(1:nbTerm_t)) = term_index;
jac_valsInit(3*nb+(1:nbTerm_t)) = 1;

jac_xInit(3*nb+nbTerm_t+(1:nbTerm_t)) = nb+term_index;
jac_yInit(3*nb+nbTerm_t+(1:nbTerm_t)) = nb+term_index;
jac_valsInit(3*nb+nbTerm_t+(1:nbTerm_t)) = 1;

currJacInd_init = 3*nb+2*nbTerm_t;

clear nbTerm_t;

% Pif jacobian
denomMut_class = zeros(1,nb);
for currClass_succ = 1:maxSucc
    % only deal with non terminal queues
    currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ))';
    % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
    if ~isempty(currQs_1)
        
        % get the indexes of their successors
        currSucc_ind = succAll(currQs_1,1:currClass_succ);
        %disp('all of the above values should be calculated during the preprocessing phase');
        
        pij_pjfull = zeros(length(currQs_1),currClass_succ);
        
        currLeng = length(currQs_1);
        for currS_id = 1:currClass_succ
            % calc the linear indexes in order to access the element:
            % pij(currQs_1,currSucc_ind(:,currS_id))
            lin_indx = (currSucc_ind(:,currS_id)-1)*nb+currQs_1;
            
            %  jacobian : Pif
            currAllVals = currJacInd_init+(1:currLeng);
            jac_xInit(currAllVals)=(nb*4+currQs_1);
            jac_yInit(currAllVals)=nb*3+currSucc_ind(:,currS_id);
            jac_valsInit(currAllVals) = -pij(lin_indx);
            currJacInd_init = currJacInd_init+currLeng;
        end;
        clear currQs_1 currSucc_ind denomMut_class1 pij_pjfull
    end;
end;


% check that nbJacInit is the exact value :
if currJacInd_init ~= nbJacInit error('recalculate the nb of jacobian elements initialized'); end;


options=optimset('Algorithm',myOptions{2},'DerivativeCheck',myOptions{4},'Display','iter','Jacobian','on','MaxFunEvals',100000000,'MaxIter',myOptions{3},'TolFun',myOptions{1},'TolX',myOptions{1}); %,'NonlEqnAlgorithm',medAlgo{1});
%,'LargeScale','on'
%'DerivativeCheck','on',


%keyboard; %save('tmpX_init_fsolve.mat','x_init');

[x_out,fval,exitflag,output,jacobian] = fsolve(@systEqtn, x_init, options);

print_fsolve_out(fval, exitflag, output,jacobian, outFile);

extra_vbles_out = zeros(nb_extra_unknowns,nb); %order mutAvg, muHat, lambda_IN, pjfull, Pif, EBi, lamOut, mut0
for j=1:nb_extra_unknowns
    extra_vbles_out(j,1:nb)=x_out((nb*(j-1)+1):(nb*j));
end;
disp('x_out contains lamEff!! not lamPiQ');
% carefull output muHat and not 1/muHat



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


    function [eq_out, jacob] = systEqtn(s0)
        %%%% coding of s0:
        % mu_t0 = s(nb+1,1:nb); unblocking rate of first blkd job.
        % mu_hat = s(nb+2,1:nb);
        % lamEFF = s(nb+3,1:nb); %INSTEAD OF lambdaPiQ
        % pj_full = s(nb+4,1:nb);
        % pj_full_succ = s(nb+5,1:nb);
        % -- P(N>=1) = P(A=1)+P(B=1) -- overdetermined syst?
        %%%%
        
        eq_out=zeros(1,nbVbles);
        % jacob initialization is done at the end (improves performance)
        jac_x = zeros(1,nbNonNullJac); % memory allocation = zeros(1,nzmax)
        jac_y = zeros(1,nbNonNullJac);
        jac_vals = zeros(1,nbNonNullJac);
        currJacInd = 0;
        
        %%%%%
        % updating of syst vbles
        %%%%
        
        
        if ~isempty(term_index) % special cases of terminal stns
            % Pif
            eq_out(nb*4+term_index)=s0(nb*4+term_index);
            
            % mut
            eq_out(term_index)=s0(term_index); % unblocking rate assumed null (CONVENTION, because theoretically it should be infinite)
            
            % mu_hat
            eq_out(nb+term_index)=s0(nb+term_index) - mu(term_index);
            
            % jacobians are calculated in the preprocessing
        end;
        
        %disp(' Cf values ci-dessous, they should be calculated during the preprocessing phase');
        % mut, Pif
        % calculate mut equations vectorized depending on the nb of successors of each queue
        denomMut_class = zeros(1,nb);
        for currClass_succ = 1:maxSucc
            % only deal with non terminal queues
            currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ))';
            % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
            if ~isempty(currQs_1)
                
                % get the indexes of their successors
                currSucc_ind = succAll(currQs_1,1:currClass_succ);
                %disp('all of the above values should be calculated during the preprocessing phase');
                denomMut_class1 = s0(nb*2+currSucc_ind)./s0(nb+currSucc_ind);
                
                % sum for a given queue across its successor queues
                % sum across ROWS
                if currClass_succ==1
                    denomMut_class(currQs_1) = denomMut_class1;
                else
                    denomMut_class(currQs_1) = sum(denomMut_class1,2);
                end;
                % full mut equation is calculated outside of this loop
                
                pij_pjfull = zeros(length(currQs_1),currClass_succ);
                
                currLeng = length(currQs_1);
                for currS_id = 1:currClass_succ
                    
                    % jacobian : mut
                    currAllVals = currJacInd+(1:currLeng);
                    jac_x(currAllVals)=currQs_1;
                    jac_y(currAllVals)=nb*2+currSucc_ind(:,currS_id);
                    jac_vals(currAllVals) = -1./s0(nb+currSucc_ind(:,currS_id));
                    currJacInd = currJacInd+currLeng;
                    
                    currAllVals = currJacInd+(1:currLeng);
                    jac_x(currAllVals)=currQs_1;
                    jac_y(currAllVals)=nb+currSucc_ind(:,currS_id);
                    jac_vals(currAllVals) = s0(nb*2+currSucc_ind(:,currS_id))./(s0(nb+currSucc_ind(:,currS_id)).^2);
                    currJacInd = currJacInd+currLeng;
                    
                    % used for Pif
                    % calc the linear indexes in order to access the element:
                    % pij(currQs_1,currSucc_ind(:,currS_id))
                    lin_indx = (currSucc_ind(:,currS_id)-1)*nb+currQs_1;
                    pij_pjfull(:,currS_id) = pij(lin_indx)'.*s0(nb*3+currSucc_ind(:,currS_id));
                    
                    %  Pif jacobian : calculated in preprocessing
                end;
                
                % Pif
                eq_out(nb*4+currQs_1) = s0(nb*4+currQs_1) - sum(pij_pjfull,2)';
                clear currQs_1 currSucc_ind denomMut_class1 pij_pjfull
            end;
        end;
        
        % mut eqtn
        eq_out(non_term_index)=s0(nb*2+non_term_index)./s0(non_term_index)-denomMut_class(non_term_index);
        
        % mut jacobian
        currLeng = length(non_term_index);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = non_term_index;
        jac_y(currAllVals)= nb*2+non_term_index;
        jac_vals(currAllVals)= 1./s0(non_term_index);
        currJacInd = currJacInd+currLeng;
        
        currLeng = length(non_term_index);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = non_term_index;
        jac_y(currAllVals) = non_term_index;
        jac_vals(currAllVals) = -s0(nb*2+non_term_index)./(s0(non_term_index).^2);
        currJacInd = currJacInd+currLeng;
        clear denomMut_class
        
        
        % muHat
        % eqtn is vectorized for single server stations
        
        denomMuHat = (1./mu(non_term_index)+(s0(nb*4+non_term_index)./s0(non_term_index)));
        eq_out(nb+non_term_index)=1./s0(nb+non_term_index)-denomMuHat;
        
        currLeng = length(non_term_index);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb+non_term_index;
        jac_y(currAllVals) = nb+non_term_index;
        jac_vals(currAllVals) = -1./(s0(nb+non_term_index).^2);
        currJacInd = currJacInd+currLeng;
        
        currLeng = length(non_term_index);
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb+non_term_index;
        jac_y(currAllVals) = non_term_index;
        jac_vals(currAllVals) = (s0(nb*4+non_term_index)./(s0(non_term_index).^2));
        currJacInd = currJacInd+currLeng;
        
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb+non_term_index;
        jac_y(currAllVals) = nb*4+non_term_index;
        jac_vals(currAllVals) = -1./s0(non_term_index);
        currJacInd = currJacInd+currLeng;
        
        
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
            currQs_p = nonSourceQ(find(nbPred(nonSourceQ)==currClass_pred))';
            % contains the q IDs of the non_term queues with nb of predecessors equal to currClass_pred
            
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
            jac_vals(currAllVals) = gamma(currQs_p);
            currJacInd = currJacInd+currLeng;
            
            clear currQs_p currPred_ind pij_pred currP lin_indx
        end;
        
        clear nonSourceQ
        
        
        % pjfull using closed form eqtns for loss models
        rho = s0(nb*2+(1:nb))./s0(nb+(1:nb));
        
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
        jac_vals(currAllVals) = -derivRho./s0(nb+(1:nb));
        currJacInd = currJacInd+currLeng;
        % wrt muHat
        currAllVals = currJacInd+(1:currLeng);
        jac_x(currAllVals) = nb*3+(1:nb);
        jac_y(currAllVals) = nb+(1:nb);
        jac_vals(currAllVals) = (derivRho.*s0(nb*2+(1:nb)))./(s0(nb+(1:nb)).^2);
        currJacInd = currJacInd+currLeng;
        clear derivRho;
        
        
        % eqtn {P(N>=1)=P(A=1)+P(B=1)}
        %disp('I believe P(N>=1) eqtn makes the system of eqtns overdetered, thus im not including it for now');
        %eq_out(nb*4+(1:nb))=  1- (1-rho)./(1-(rho.^(buffer+2))) - ( ( s0(nb*2+(1:nb))./(s0(nb*2+(1:nb))+s0(nb+(1:nb))) ).*( 1+ (allMu.*s0(nb*4+(1:nb)))./(s0(1:nb)+allMu.*s0(nb*4+(1:nb)).*(s0(nb*2+(1:nb))./(s0(nb*2+(1:nb))+s0(nb+(1:nb)))) ) );
        %afaf derivs
        
        
        % check that nbJacInit is the exact value :
        if currJacInd ~= length(jac_x) error('recalculate the nb of jacobian elements initialized or copy jac_x (code follows)'); end;
        
        jacob = sparse([jac_x,jac_xInit],[jac_y,jac_yInit],[jac_vals,jac_valsInit],nbVbles,nbVbles);
        % used for fmincon
        % UNLIKE fsolve:  fmincon defines jacob(i,j) is the partial derivative of equation(j) with respect to x(i)
        % thus need to transpose jacob_2dim;
        %jacob_2dim=jacob_2dim';
        
        %grad_ceq=[jacob_2dim;zeros(1,nbVbles)]; % adding derivatives wrt to z (they are all null)
        %max(abs(eq_out))
        %cond(jacob)
        %condest(jacob) % for sparse structure
        %keyboard; % save('simLo_fsolve.mat');
        
        
        clear jac_vals jac_x jac_y s0 currAllVals %j_tmp; %jacob_2dim;
    end % fctn nonLinConstr
end % fctn netwk1eqtn