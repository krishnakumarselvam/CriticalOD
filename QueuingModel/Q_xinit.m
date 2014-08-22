%Modified ny Krishna on Nov 21 2013
function xOut = Q_xinit(scenarios, serv,pij,mu,lambda_flowCvtn,nb,nb_extra_unknowns,term_index,non_term_index,maxSucc,nbSucc,succAll,buffer,gamma)


xOut = zeros(1,nb*nb_extra_unknowns);
Rho_hat_index = 1:nb;
Rho_eff_index = nb+1:2*nb;
Lambda_eff_index = 2*nb+1:3*nb;
Pnk_index=3*nb+1:4*nb;
Pi_full_index = 4*nb+1:5*nb;
true_rho_index=5*nb+1:6*nb;


currCase = 0;

%Reducing the lambda in over-congested links
%Finding the probability of Ni = Ki , using a uniform approximation.

rho=lambda_flowCvtn./mu;
%rho=rho/(max(rho))*0.9;
for init_iter = 1:1
xOut(Pnk_index)=((1-rho).*(rho.^(buffer)))./(1-(rho.^(buffer+1)));
xOut(xOut==1)=0.99;
temp = xOut(Pnk_index);
AA = eye(nb)-pij';

bb=gamma.*(1-xOut(Pnk_index));
right_lambda = AA\bb';
xOut(Lambda_eff_index)=right_lambda;
disp(sum(xOut(Lambda_eff_index)))


% figure
% bar(bb);
% figure
% bar(right_lambda,'r')

%compute Pi full, no need to change
if currCase == 0
    %Pi full for terminal queues
    xOut(nb*4+term_index) = 0;
    % non_terminal
    for currClass_succ = 1:maxSucc
        % only deal with non terminal queues
        currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ))';
        % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
        if ~isempty(currQs_1)
            % get the indexes of their successors
            currSucc_ind = succAll(currQs_1,1:currClass_succ);
            pij_pjfull = zeros(length(currQs_1),currClass_succ);
            for currS_id = 1:currClass_succ
                % calc the linear indexes in order to access the element:
                % pij(currQs_1,currSucc_ind(:,currS_id))
                lin_indx = (currSucc_ind(:,currS_id)-1)*nb+currQs_1;
                pij_pjfull(:,currS_id) = pij(lin_indx)'.*xOut(nb*3+currSucc_ind(:,currS_id));
            end;
            xOut(nb*4+currQs_1) = sum(pij_pjfull,2)';
        end;
    end;
     
   % vbles are lambda(eff)/mut and lambda(eff)/muHat thus we simply need to solve a linear system of equations
    nbNonTerm = length(non_term_index);
    nonNullElem = sum(nbSucc)+2*nb+nbNonTerm;
    rowsA = zeros(1,nonNullElem);
    colsA = zeros(1,nonNullElem);
    valsA = zeros(1,nonNullElem);

    % diagonal elements (coeff of rho_tilda)
    rowsA(1:nb) = nb+(1:nb);
    colsA(1:nb) = nb+(1:nb);
    valsA(1:nb) = 1;
    currElem=nb;

    currAllVals = currElem+(1:length(term_index));
    rowsA(currAllVals) = term_index;
    colsA(currAllVals) = term_index;
    valsA(currAllVals) = 1;
    currElem = currElem + length(term_index);

    % diagonal elements for non_term xt
    currAllVals = currElem+(1:nbNonTerm);
    rowsA(currAllVals) = non_term_index;
    colsA(currAllVals) = non_term_index; % 
    valsA(currAllVals) = 1;%xt diagonal coeffients are 1
    currElem = currElem + nbNonTerm;

    if currElem~=(2*nb) error('check'); end;


    % Rho effective equation for nonTerm queues
    currAllVals = currElem+(1:nbNonTerm);
    rowsA(currAllVals) = nb + non_term_index;
    colsA(currAllVals) = non_term_index;
    valsA(currAllVals) = -xOut(nb*4+non_term_index);
    currElem = currElem + nbNonTerm;

    % rho tilda eqtn for nonTerm stns
    [tmpV,currQs,currSuccs]=find(succAll');
    if length(currSuccs)~=sum(nbSucc) error('incoherent total number of successors'); end;

    currAllVals = currElem+(1:length(currSuccs));
    rowsA(currAllVals) = currQs;
    colsA(currAllVals) = nb + currSuccs; % muHat of the successor stns
    valsA(currAllVals) = -1;% New mut equation 
    currElem = currElem + length(currSuccs);
    A = sparse(rowsA,colsA,valsA,2*nb,2*nb);
    b = zeros(2*nb,1);
    b(Rho_eff_index) = xOut(Lambda_eff_index)./mu;
    xs = A \ b;%%xt and xhat 
    currIndx = find(xs);%associated indices
    xOut(currIndx) = xs(currIndx);%update xout
    
    xOut(true_rho_index)=xOut(Rho_eff_index)./(1-xOut(Pnk_index));
 
    NUM = 322;
%     indices = [NUM,nb+NUM,3*nb+NUM,4*nb+NUM,5*nb+NUM];
%     bar(xOut(indices))
%     disp(xOut(nb*2+NUM))
%     close all
%     scatter(xOut(true_rho_index),rho);
%     hold on
%     plot(xOut(true_rho_index),xOut(true_rho_index),'r');
%    
end;
    if(init_iter==1)
         x_old = xOut;
         rho=x_old(true_rho_index);
    else
         x_old = init_iter/(init_iter+1)* x_old + 1/(init_iter+1)*xOut;  
         rho=x_old(true_rho_index);
    end
end

xOut=x_old;

loc=find(isnan(xOut));
end