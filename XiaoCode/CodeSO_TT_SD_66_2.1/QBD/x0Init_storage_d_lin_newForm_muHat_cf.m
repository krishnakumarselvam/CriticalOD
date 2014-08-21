function xOut = x0Init_storage_d_lin_newForm_muHat_cf(scenarios, serv,pij,mu,lambda_flowCvtn,nb,nb_extra_unknowns,term_index,non_term_index,maxSucc,nbSucc,succAll,buffer)
% differs from x0Init_storage_d_lin_newForm_muHat: full dbns are no longer caculated (closed form expression of the loss models is used for pfull)


if ~isempty(find(serv~=1)) error('use code x0Init_storage_d_lin_Form_old for  non single server networks (that code is still to be derivChecked and debugge'); end;

xOut = zeros(1,nb*nb_extra_unknowns);

currCase = 0;

% pjfull
switch scenarios{1}
    
    case 'piRdm'
        rand('state', sum(100*clock)); % Initialize rand to a different state each time.
        xOut(nb*3+(1:nb)) = rand(1,nb);
        
    case 'piUnif'
        xOut(nb*3+term_index) = 1./(buffer(term_index)+2);
        xOut(nb*3+non_term_index) = 2./(2.*buffer(non_term_index)+3);
        
    case 'piCFApprox' % closed form expression using mu instead of muHat
        rhoTmp = lambda_flowCvtn./mu;
        xOut(nb*3+(1:nb)) = ((1-rhoTmp).*(rhoTmp.^(buffer+1)))./(1-(rhoTmp.^(buffer+2)));
        clear rhoTmp;
        
    case 'piLamMuHat' %fix lambda and muHat, deduce all other parameters
        % do nothing for now, but avoid error
        currCase = 1;
        
    otherwise error('pi initialization not identified');
end;

if currCase == 0
    
    % lambda
    xOut(nb*2+(1:nb)) = lambda_flowCvtn;
    
    % Pif
    % terminal
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
    
    
    % mu
    switch scenarios{2}
        case 'Mu'
            xOut(1:nb) = mu;
            xOut(nb+(1:nb)) = mu;
            
        case 'DeducedEqtn'
            
            % vbles are 1/mut and 1/muHat thus we simply need to solve a linear system of equations
            nbNonTerm = length(non_term_index);
            nonNullElem = sum(nbSucc)+2*nb+nbNonTerm;
            rowsA = zeros(1,nonNullElem);
            colsA = zeros(1,nonNullElem);
            valsA = zeros(1,nonNullElem);
            
            % diagonal elements for all muHat and terminal mut
            rowsA(1:nb) = nb+(1:nb);
            colsA(1:nb) = nb+(1:nb);
            valsA(1:nb) = 1;
            currElem=nb;
            
            currAllVals = currElem+(1:length(term_index));
            rowsA(currAllVals) = term_index;
            colsA(currAllVals) = term_index;
            valsA(currAllVals) = 1;
            currElem = currElem + length(term_index);
            
            % diagonal elements for non_term mut
            currAllVals = currElem+(1:nbNonTerm);
            rowsA(currAllVals) = non_term_index;
            colsA(currAllVals) = non_term_index; % muHat of the successor stns
            valsA(currAllVals) = (xOut(nb*2+non_term_index).*(1-xOut(nb*3+non_term_index)));
            currElem = currElem + nbNonTerm;
            
            if currElem~=(2*nb) error('check'); end;
            
            
            % muHat eqtn for nonTerm stns
            currAllVals = currElem+(1:nbNonTerm);
            rowsA(currAllVals) = nb + non_term_index;
            colsA(currAllVals) = non_term_index;
            valsA(currAllVals) = -xOut(nb*4+non_term_index);
            currElem = currElem + nbNonTerm;
            
            % mut eqtn for nonTerm stns
            [tmpV,currQs,currSuccs]=find(succAll');
            if length(currSuccs)~=sum(nbSucc) error('incoherent total number of successors'); end;
            
            currAllVals = currElem+(1:length(currSuccs));
            rowsA(currAllVals) = currQs;
            colsA(currAllVals) = nb + currSuccs; % muHat of the successor stns
            valsA(currAllVals) = -(xOut(nb*2+currSuccs).*(1-xOut(nb*3+currSuccs)));
            currElem = currElem + length(currSuccs);
            
            % OLD vsn BUG: valsA(currAllVals) = -(xOut(first_index+nb*2+currSuccs-1).*(1-xOut(first_index+nb*3+currSuccs-1)))./(xOut(first_index+nb*2+currQs-1).*(1-xOut(first_index+nb*2+currQs-1))); % because serv==1
            
            A = sparse(rowsA,colsA,valsA,2*nb,2*nb);
            
            b = zeros(2*nb,1);
            b(nb+(1:nb)) = 1./mu;
            
            musInv = A \ b;
            currIndx = find(musInv);
            xOut(currIndx) = 1./musInv(currIndx);
            
            return;
            
        otherwise error('mu initialization not identified');
    end;
    
else % currCase == 1 : fix lambda and muHat, deduce all other parameters
    % check
    if (currCase~=1) | (~strcmp(scenarios{2},'muLamMuHat')) | (~strcmp(scenarios{1},'piLamMuHat'))
        error('case not defined');
    end;
    
    % lambda, muHat
    xOut(nb*2+(1:nb)) = lambda_flowCvtn;
    xOut(nb+(1:nb)) = mu;
    
    %pifull
    rhoTmp = xOut(nb*2+(1:nb))./xOut(nb+(1:nb));
    xOut(nb*3+(1:nb)) = ((1-rhoTmp).*(rhoTmp.^(buffer+1)))./(1-(rhoTmp.^(buffer+2)));
    
    % mut, Pif
    denomMut_class = zeros(1,nb);
    for currClass_succ = 1:maxSucc
        % only deal with non terminal queues
        currQs_1 = non_term_index(find(nbSucc(non_term_index)==currClass_succ))';
        % contains the q IDs of the non_term queues with nb of successors equal to currClass_succ
        if ~isempty(currQs_1)
            
            % get the indexes of their successors
            currSucc_ind = succAll(currQs_1,1:currClass_succ);
            %disp('all of the above values should be calculated during the preprocessing phase');
            denomMut_class1 = (xOut(nb*2+currSucc_ind).*(1-xOut(nb*3+currSucc_ind)))./xOut(nb+currSucc_ind);
            
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
                
                % used for Pif
                % calc the linear indexes in order to access the element:
                % pij(currQs_1,currSucc_ind(:,currS_id))
                lin_indx = (currSucc_ind(:,currS_id)-1)*nb+currQs_1;
                pij_pjfull(:,currS_id) = pij(lin_indx)'.*xOut(nb*3+currSucc_ind(:,currS_id));
                
                %  Pif jacobian : calculated in preprocessing
            end;
            
            % Pif
            xOut(nb*4+currQs_1) = sum(pij_pjfull,2)';
            clear currQs_1 currSucc_ind denomMut_class1 pij_pjfull
        end;
    end;
    
    % mut eqtn
    mutInv = denomMut_class(non_term_index)./(xOut(nb*2+non_term_index).*(1-xOut(nb*3+non_term_index)));
    xOut(non_term_index) = 1./mutInv;
end;