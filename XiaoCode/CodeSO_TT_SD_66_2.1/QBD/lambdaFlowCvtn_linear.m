function lambda_flowCvtn = lambdaFlowCvtn_linear(gamma, pij,outFile)
% differs from lambdaFlowCvtn in that eqtns are defined as a linear system + solved with / operator instead of fsolve
% calculates the flow conservation arrival rates


% INPUTS
% gamma             vector of external arrival rates (of length number of stations)
% pij               matrix of routing/transition probabilities
% OUTPUT
% lambda_flowCvtn: row vector

[nb,nb_tmp] = size(pij);

[r,c]=size(gamma);
if c~=1 % gamma is a row vector
    gamma_local = gamma';
else
    gamma_local = gamma;
end;

eye_sparse = sparse(1:nb,1:nb,1,nb,nb);
lambda_flowCvtn = ((eye_sparse-pij')\gamma_local)';

% printing output
if ~isempty(outFile)
    fid1 = fopen(outFile,'a');
    %keyboard;
    fprintf(fid1,'%s\n','lambda_flow_cvtn equals');
    for i = 1:nb
        fprintf(fid1,'%1.6f\n',lambda_flowCvtn(i));
    end;
    fprintf(fid1,'\n\n');
    fclose(fid1);
end;
