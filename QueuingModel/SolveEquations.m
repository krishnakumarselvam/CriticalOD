
function [x_out,fval,exitflag,output,jacobian,extra_vbles_out] = SolveEquations(gamma, lambda_flowCvtn, mu, serv, buffer, pij,scenarios,myOptions,knownX0,loadX0,iter_rdm)
extra_vbles_out=0;
nb = length(buffer); %nb of nodes
Rho_hat_index = 1:nb;
Rho_eff_index = nb+1:2*nb;
Lambda_eff_index = 2*nb+1:3*nb;
Pnk_index=3*nb+1:4*nb;
Pi_full_index = 4*nb+1:5*nb;
true_rho_index=5*nb+1:6*nb;

[term_index,non_term_index,nbSucc,maxSucc,succAll,nbPred,maxPred,predAll] = prelim_storage_cf(serv, pij, buffer);

nb_extra_unknowns = 6; % does not include mu nor phase vbles.
nbVbles = nb*nb_extra_unknowns;

%Defining Normalizing factors

NORMALIZING_FACTOR = ones(1,nbVbles);
NORMALIZING_FACTOR(Lambda_eff_index)=1;

disp('calculating x0 for Queuing Model');
% in this initialization mu is considered exogenous. i.e we initialize vbles to satisfy state eqtns but decision vbles play no role in the initialization.
x_init = Q_xinit(scenarios, serv,pij,mu,lambda_flowCvtn,nb,nb_extra_unknowns,term_index,non_term_index,maxSucc,nbSucc,succAll,buffer,gamma);
disp('x0 calculation done for Queuing Model');


sumSucc = sum(nbSucc);
sumPred = sum(nbPred);
%nbNonSourceQ = find(nbPred>0);
nbNonNullJac = nb + length(non_term_index)*2  + length(find(nbPred>0))*1+length(find(nbPred==0));
%nbNonNullJac
clear sumSucc sumPred;


%%%%%%
% calculate all the jacobian values that are constant.
%%%%%%


Jac_X0=1:nbVbles;
Jac_Y0=1:nbVbles;
Jac_Val0=ones(1,nbVbles);


%Equation 1 - constant jacobian terms
Jac_X1=zeros(0);
Jac_Y1=zeros(0);
for currClass_succ = 1:maxSucc
    % only deal with non terminal queues
    currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ))';
    % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
    if ~isempty(currQs_1)
        % get the indexes of their successors
        currSucc_ind = succAll(currQs_1,1:currClass_succ);
        
        for currS_id = 1:currClass_succ
          l=length(Jac_X1);
          Jac_X1(l+(1:length(currQs_1)))= currQs_1';
          Jac_Y1(l+(1:length(currQs_1)))=nb + currSucc_ind(:,currS_id);
          Jac_Val1(l+(1:length(currQs_1)))=-1;
        end;
        
    end;
end;

%Equation 2 - constantjacobian terms
Jac_X2_init=Rho_eff_index; 
Jac_Y2_init=Lambda_eff_index; 
Jac_Val2_init=(-1)./mu;


nbTerm_t = length(term_index);
nbnonTerm_t=length(non_term_index);
nbJacInit = 3*nb+3*nb+2*sum(nbSucc);
jac_xInit = zeros(1,nbJacInit); % memory allocation = zeros(1,nzmax)
jac_yInit = zeros(1,nbJacInit);
jac_valsInit = zeros(1,nbJacInit);

% diagonal elements lambda, pjfull, Pif
jac_xInit(1:(3*nb)) = (nb*2) + (1:(3*nb));
jac_yInit(1:(3*nb)) = (nb*2) + (1:(3*nb));
jac_valsInit(1:(3*nb)) = 1;
% diagonal elements for term Qs: rhot, rhoHat
jac_xInit(3*nb+(1:nbTerm_t)) = term_index;
jac_yInit(3*nb+(1:nbTerm_t)) = term_index;
jac_valsInit(3*nb+(1:nbTerm_t)) = 1;

jac_xInit(3*nb+nbTerm_t+(1:nbTerm_t)) = nb+term_index;
jac_yInit(3*nb+nbTerm_t+(1:nbTerm_t)) = nb+term_index;
jac_valsInit(3*nb+nbTerm_t+(1:nbTerm_t)) = 1;

currLength=3*nb+2*nbTerm_t;
%xt jacobians diagonal term
jac_xInit(currLength+(1:nbnonTerm_t)) = non_term_index;
jac_yInit(currLength+(1:nbnonTerm_t)) = non_term_index;
jac_valsInit(currLength+(1:nbnonTerm_t)) = 1;

% xhat jacobians diagonal term
currLength=3*nb+2*nbTerm_t+nbnonTerm_t;
jac_xInit(currLength+(1:nbnonTerm_t)) = nb+non_term_index;
jac_yInit(currLength+(1:nbnonTerm_t)) = nb+non_term_index;
jac_valsInit(currLength+(1:nbnonTerm_t)) = 1;

currJacInd_init = 3*nb+2*nb;


clear nbTerm_t;

%  xt jacobian 
% xt_class = zeros(1,nb);
for currClass_succ = 1:maxSucc
    % only deal with non terminal queues
    currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ))';
    % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
    if ~isempty(currQs_1)
        
        % get the indexes of their successors
        currSucc_ind = succAll(currQs_1,1:currClass_succ);
        
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
            
            %jacobian :xt
            currAllVals = currJacInd_init+(1:currLeng);
            jac_xInit(currAllVals)=currQs_1;
            jac_yInit(currAllVals)=nb+currSucc_ind(:,currS_id);
            jac_valsInit(currAllVals) = -1;
            currJacInd_init = currJacInd_init+currLeng;
        end;
        clear currQs_1 currSucc_ind denomMut_class1 pij_pjfull
    end;
end;

%rho eff Jacobian (wrt lambda effective)
currAllVals = currJacInd_init+(1:nb);
jac_xInit(currAllVals)=nb+(1:nb);
jac_yInit(currAllVals)=nb*2+(1:nb);
jac_valsInit(currAllVals) = -1./mu;
 
%Linsen: not sure about is currently, need to check
currJacInd_init = currJacInd_init+nb;
% check that nbJacInit is the exact value :
if currJacInd_init ~= nbJacInit error('recalculate the nb of jacobian elements initialized'); end;
%'MaxIter',myOptions{3}
%options=optimset('Algorithm',myOptions{2},'DerivativeCheck',myOptions{4},'Display','on','Jacobian','on','MaxFunEvals',100000000,'MaxIter',myOptions{3},'TolFun',myOptions{1},'TolX',myOptions{1}); %,'NonlEqnAlgorithm',medAlgo{1});
options=optimset('Algorithm',myOptions{2},'DerivativeCheck','off','Display','off','Jacobian','on','MaxFunEvals',100000000,'MaxIter',100,'TolFun',myOptions{1},'TolX',myOptions{1}); 

[EQ,Jacobian]=systEqtn(x_init);

[x_out,fval,exitflag,output,jacobian] = fsolve(@systEqtn, x_init, options);
close all


    function [eq_out, jacob] = systEqtn(s0)
       eq_out=zeros(1,nbVbles);
        % jacob initialization is done at the end (improves performance)
        %jacob=zeros(nbVbles,nbVbles);
        
                
        
        if ~isempty(term_index) % special cases of terminal queues
            % Pif, same
            eq_out(nb*4+term_index)=s0(nb*4+term_index);
            % xt, Linsen: I assume xt=zero for non terminal?
            eq_out(term_index)=s0(term_index); % unblocking rate assumed null (CONVENTION, because theoretically it should be infinite)
            % x_hat
            eq_out(nb+term_index)=s0(nb+term_index) - s0(2*nb+term_index)./mu(term_index);
        end;
        
        %Equation 5 and a bit of Equation 1
        

        Jac_X5_part2=zeros(0);
        Jac_Y5_part2=zeros(0);
        Jac_Val5_part2=zeros(0);
         l=0;
        xt_class = zeros(1,nb);
        for currClass_succ = 1:maxSucc
            % only deal with non terminal queues
            currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ))';
            % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
            if ~isempty(currQs_1)
                % get the indices of their successors
                currSucc_ind = succAll(currQs_1,1:currClass_succ);
                % currSucc_ind is a 2D array
                xt_class1 = s0(nb+currSucc_ind);
                % sum for a given queue across its successor queues
                % sum across ROWS
              
                if currClass_succ==1
                    xt_class(currQs_1) = xt_class1;
                else
                    xt_class(currQs_1) = sum(xt_class1,2);
                end;
                % full mut equation is calculated outside of this loop
                pij_pjfull = zeros(length(currQs_1),currClass_succ);
                
                currLeng = length(currQs_1);
               
                for currS_id = 1:currClass_succ
                    lin_indx = (currSucc_ind(:,currS_id)-1)*nb+currQs_1;
                    pij_pjfull(:,currS_id) = pij(lin_indx)'.*s0(nb*3+currSucc_ind(:,currS_id));
                    
                     
                    Jac_X5_part2(l+(1:currLeng))=nb*4+currQs_1;
                    Jac_Y5_part2(l+(1:currLeng))=nb*3+currSucc_ind(:,currS_id);
                    Jac_Val5_part2(l+(1:currLeng))=-pij(lin_indx);
                    l=l+currLeng;
                    
                end;
                
                
                
                % Pif equation No.5
                eq_out(nb*4+currQs_1) = s0(nb*4+currQs_1) - sum(pij_pjfull,2)';%Linsen: no change on Equation 5
                clear currQs_1 currSucc_ind tt_class1 pij_pjfull
            end;
        end;
        
        Jac_X5=[ Jac_X5_part2];
        Jac_Y5=[ Jac_Y5_part2];
        Jac_Val5=[ Jac_Val5_part2];
        
        
        
        %Equation 1: Rho tilda 
        eq_out(non_term_index)=s0(non_term_index)-xt_class(non_term_index);%Linsen Change
        
        
        %Equation 2: Rho effective
        eq_out(nb+non_term_index)=s0(nb+non_term_index)-s0(nb*2+non_term_index)./mu(non_term_index)-s0(nb*4+non_term_index).*s0(non_term_index);
        
        %Equation 2 Jacobian (refer to the beginning of file for constant
        %terms
        Jac_X2=[Rho_eff_index,Rho_eff_index];
        Jac_Y2=[Pi_full_index, Rho_hat_index];
        Jac_Val2=[-s0(Rho_hat_index),-s0(Pi_full_index)];
        
        %Equation 3: Lambda effective
        sourceQ = find(nbPred==0);
        eq_out(nb*2+sourceQ)=s0(nb*2+sourceQ) - gamma(sourceQ).*(1-s0(nb*3+sourceQ));
       
        %Jacobian for source Queues
        Jac_X3_src = [ Lambda_eff_index(sourceQ)];
        Jac_Y3_src = [Pnk_index(sourceQ) ];
        Jac_Val3_src = [(gamma(sourceQ))];
        
        nonSourceQ = find(nbPred>0);
        

        Jac_X3_nonsrc_part3=zeros(0);
        Jac_Y3_nonsrc_part3=zeros(0);
        Jac_Val3_nonsrc_part3=zeros(0);
        
        for currClass_pred = 1:maxPred
            
            currQs_p = nonSourceQ(find(nbPred(nonSourceQ)==currClass_pred))';
            % contains the q IDs of the non_term queues with nb of predecessors equal to currClass_pred
            % get the indexes of their predecessors
            if(length(currQs_p))
            currPred_ind = predAll(currQs_p,1:currClass_pred);
            num_q=length(currQs_p);
            pij_pred = zeros(length(currQs_p),currClass_pred);
            
            for currP = 1:currClass_pred
                % get the linear index of the elements
                % pij(currPred_ind(:,currP),currQs_p)
                lin_indx = ((currQs_p - 1)*nb + currPred_ind(:,currP))';
                pij_pred(:,currP) = pij(lin_indx).*s0(nb*2+currPred_ind(:,currP));
                l=length(Jac_X3_nonsrc_part3);
                Jac_X3_nonsrc_part3(l+(1:num_q))=currQs_p + 2*nb;
                Jac_Y3_nonsrc_part3(l+(1:num_q))=2*nb + currPred_ind(:,currP);
                Jac_Val3_nonsrc_part3(l+(1:num_q))=-pij(lin_indx);
                
            end;
            % lambda eqtn.
            currSum = (sum(pij_pred,2))';
            eq_out(nb*2+currQs_p) = s0(nb*2+currQs_p) - (gamma(currQs_p).*(1-s0(nb*3+currQs_p)) + currSum);
            end
        end;
        Jac_X3_nonsrc_part1 = [2*nb + nonSourceQ];
        Jac_Y3_nonsrc_part1 = [Pnk_index(nonSourceQ)];
        Jac_Val3_nonsrc_part1 = [(gamma(nonSourceQ))];
        
        Jac_X3 = [Jac_X3_src, Jac_X3_nonsrc_part1,Jac_X3_nonsrc_part3];
        Jac_Y3 = [Jac_Y3_src, Jac_Y3_nonsrc_part1,Jac_Y3_nonsrc_part3];
        Jac_Val3 = [Jac_Val3_src, Jac_Val3_nonsrc_part1,Jac_Val3_nonsrc_part3];
    
        clear nonSourceQ
        
        %Equation 4 : P(Ni = Ki)
        rho = s0(true_rho_index);
        eq_out(nb*3+(1:nb))= s0(nb*3+(1:nb)) - ((1-rho).*(rho.^(buffer)))./(1-(rho.^(buffer+1)));
       
        
        
        %Jacobian for equation 4
        %derivRho = ((rho.^(buffer-1)).*((buffer-1)+1-((buffer-1)+2).*rho))./(1-(rho.^((buffer-1)+2)))+(((buffer-1)+2).*(rho.^(2.*((buffer-1)+1))).*(1-rho))./((1-(rho.^((buffer-1)+2))).^2); % Check this out...
        new_deriv = rho.^buffer./(rho.^(buffer + 1) - 1)        +      (buffer.*rho.^(buffer - 1).*(rho - 1))./(rho.^(buffer + 1) - 1)       -      (rho.^(2*buffer).*(buffer + 1).*(rho - 1))./(rho.^(buffer + 1) - 1).^2;
        
        Jac_X4=nb*3 + (1:nb);
        Jac_Y4=true_rho_index;
        Jac_Val4= -new_deriv;
        
        eq_out(true_rho_index)=s0(true_rho_index)-s0(Rho_eff_index)./(1-s0(Pnk_index));
        Jac_X14=[true_rho_index ,true_rho_index];
        Jac_Y14=[Rho_eff_index , Pnk_index];
        
        Jac_Val14=[ -1./(1-s0(Pnk_index)) , -s0(Rho_eff_index)./(1-s0(Pnk_index)).^2];
        
        JAC_X=[Jac_X0,Jac_X1,Jac_X2_init,Jac_X2,Jac_X3,Jac_X4,Jac_X5,Jac_X14,jac_xInit];
        JAC_Y=[Jac_Y0,Jac_Y1,Jac_Y2_init,Jac_Y2,Jac_Y3,Jac_Y4,Jac_Y5,Jac_Y14,jac_yInit];
        JAC_VAL=[Jac_Val0,Jac_Val1,Jac_Val2_init,Jac_Val2,Jac_Val3,Jac_Val4,Jac_Val5,Jac_Val14,jac_valsInit];
        
        jacob = sparse(JAC_X,JAC_Y,JAC_VAL,nbVbles,nbVbles);
%         for iii=1:nbVbles
%             jacob(:,iii)=jacob(:,iii)./NORMALIZING_FACTOR';
%         end
%      
%         eq_out = eq_out./NORMALIZING_FACTOR;
        
%        if(DISP_BAR)
%             close all
%             figure
%             bar(eq_out);
%             hold on
%             bar(s0,'r')
% %             figure
% %             imagesc(jacob)
% %             temp=full(jacob);
% %             colorbar
%             temp = jacob(true_rho_index,:);
%             keyboard
%        end

    end % fctn nonLinConstr

end % fctn netwk1eqtn