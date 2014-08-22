function [term_index,non_term_index,nbSucc,maxSucc,succAll,nbPred,maxPred,predAll] = prelim_storage_cf(serv, pij, buffer)
% differs from prelim_storage_newForm: block infos are no longer calculates (no longer needed since the dbns are no longer caculated (closed form expression of the loss models is used for pfull))

if ~isempty(find(serv>1)) error('this code is for single server networks'); end;
clear serv;
% INPUTS
% serv      vector of number of servers per station
% buffer    vector of buffer sizes per stn
% pij       matrix of routing probabilities

%%%%
% succ: matrix containing the succeeding nodes of node j in column j.
nb=length(buffer);
term_index = [];
non_term_index = [];
succ = {};
nbSucc = zeros(1,nb);
maxSuc_guess = 9;
succAll = zeros(nb,maxSuc_guess); % this storage could be improved since the nb of queues with nbSucc>2 is small
%succAll : matrix of size (nb,maxNbSucc )
%such that element succAll(q,1:nbSucc(q)) has the indexes of the successors of queue q.
pred = {};
nbPred = zeros(1,nb);
maxPred_guess = 11;
predAll = zeros(nb,maxPred_guess); % this storage could be improved since the nb of queues with nbSucc>2 is small

for stn =1:nb
    succ{stn}=find(pij(stn,:)~=0);
    if isempty(succ{stn}) %stn is terminal
        term_index=[term_index,stn];
    else
        non_term_index=[non_term_index,stn];
    end;
    nbSucc(stn) = length(succ{stn});
    succAll(stn,1:length(succ{stn})) = succ{stn};

    % predecessors of stn
    pred{stn}=find(pij(:,stn)~=0);
    nbPred(stn) = length(pred{stn});
    predAll(stn,1:length(pred{stn})) = pred{stn};
end

clear succ pred;

% resize succAll if maxSuc_guess was overestimated
maxSucc = max(nbSucc);
if maxSucc < maxSuc_guess
    succAll(:,(maxSucc+1):maxSuc_guess) =[];
else
    if maxSucc > maxSuc_guess
        error(' we have left out some successors');
    end;
end;

% idem for maxPred_guess
maxPred = max(nbPred);
if maxPred < maxPred_guess
    predAll(:,(maxPred+1):maxPred_guess) =[];
else
    if maxPred > maxPred_guess
        error(' we have left out some predecessors');
    end;
end;
