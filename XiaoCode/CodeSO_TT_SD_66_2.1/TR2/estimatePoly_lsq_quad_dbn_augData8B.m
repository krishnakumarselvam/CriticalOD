function betaFile_out = estimatePoly_lsq_quad_dbn_augData8B(ptDim,zFile,yFile,xkCenter_ID,currIter,allBeta_file,dbnIter_s,betaF_prev,initPtIdS,dirName,fileNwk,phases2Remove)
% differs from estimatePoly_lsq_quad_dbn_augData7: 
% 1- initial value for alpha is 1
%
% phi is now a linear poly instead of a quadratic given the observations (z, y(z)), estimate the model alphaOpt*fQ + Phi(z) by least squares minimization,
% where Phi(z) = aOpt + bOpt*x + COpt*x^2
% - z : green times
% - y(z) : average travel time
% - ptDim : dimension of decision variables (ie nb of endogenous phases)

betaF_str = [dirName,'Estimates/beta_quad_',initPtIdS,'_',dbnIter_s,'_'];

currStr = '';
for j=1:(ptDim+2)
    currStr =[currStr,'%f '];
end;


fid=fopen(zFile,'r');
% columns: ptID, pt
%dataZ_cell = textscan(fid,['%d ',currStr]);
dataZ_cell = textscan(fid,currStr);
fclose(fid);
dataZ = cell2mat(dataZ_cell);
clear dataZ_cell;

disp([int2str(length(dataZ(:,1))),' obsvtn read']);

fid=fopen(yFile,'r');
% columns: ptID replicID avgDelay
dataY = textscan(fid,'%d %d %f');
fclose(fid);

% check that (xk, f(xk)) have been stored in the observation files
if isempty(find(dataZ(:,1)==xkCenter_ID)) error('current ID has not been included in the file of observations'); end;
if isempty(find(dataY{1}==xkCenter_ID)) error('current ID has not been included in the file of observations'); end;

% check that order of observations (ie sequence of point IDs) is consistent
if (dataZ(:,1)~=dataY{1}) error('pt and observation files seem to be in a different order or of different length'); end;

[nbObs,c_tmp] = size(dataZ);
if (c_tmp ~= (ptDim+2)) error('mismatch between dimensions'); end;

currLine = find(dataZ(:,1)==xkCenter_ID);
if (dataY{1}(currLine)~=xkCenter_ID) error('ordering differs between files'); end;
xkCenter = dataZ(currLine,3:end);
% this error should never occur since the ordering was checked before
%xkAll = repmat(xk,nbObs,1);
%if length(z)~=length(xkAll) error('mismatch of the sizes of z and xkAll'); end;


% augmented data to ensure that during the initial runs matrix is of full rank
weightAugData = 10^-1; %will then be summed so will be of the order 10^-2

% CAREFULL if these indices are changed, manually adapt the jacobian indices (those within the loop)

load(fileNwk,'phase2Node');
nbNodes = length(unique(phase2Node));
aIndx = 1;
bIndx = 1+(1:(ptDim-nbNodes)); % lin terms
cIndx = bIndx(end)+(1:(ptDim-nbNodes)); 
%bIndx = 2:(ptDim+1); % lin terms
%cIndx = (ptDim+1)+(1:ptDim); 
nbVblesBeta = cIndx(end);
clear phase2Node nbNodes;

if currIter>0
    load(betaF_prev,'betaOpt');
    beta0 = betaOpt;
    clear betaOpt;
else
    beta0 = zeros(nbVblesBeta,1);
end;

options=optimset('Diagnostics','on','LargeScale','off','MaxFunEvals',100000000,'MaxIter',100000000,...
        'TolFun',1e-7,'DerivativeCheck','off'); 
% 'JacobMult'

% remove iteration indexes to facilitate manipulation of Z matrix.
dataZ(:,1:2) = [];
distAll = dataZ-repmat(xkCenter,nbObs,1);
norm_tmp = sqrt(sum((distAll').^2));
weightAll = 1./(1+norm_tmp)';
%weightAll = 1./(1+norm_tmp.^2);
%weightAll = (1./(1+sum((distAll.^2)').^2))';
clear distAll;

figure;
stem([weightAugData*ones(1,nbVblesBeta),weightAll']);
title(['weights for augmentedData + per observation. Iter:',int2str(currIter)]);


% remove one phase per intersection. 
% phases to remove were chosen based on redunduncies with other phases
C1 = zeros(nbObs,nbVblesBeta);
dataZ2 = dataZ;
clear dataZ;
dataZ2(:,phases2Remove)=[];
C1 = [ones(nbObs,1),dataZ2,dataZ2.^2];
% multiply by their weights
for j=1:nbVblesBeta
    C1(:,j) = weightAll.*C1(:,j);
end;
C = [C1;weightAugData*eye(nbVblesBeta)];
clear C1;
d = [weightAll.*dataY{3}; weightAugData*[zeros(nbVblesBeta,1)]]; %augData sets betas to zero and aplha to 1


[betaOpt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,[],[],[],[],[],[],beta0,options);

betaFile_out = [betaF_str,int2str(currIter),'.mat'];
save(betaFile_out);

% print to outFile
% identifiers might not be unique
currStr = [];
for j = 1:nbVblesBeta
    currStr = [currStr,'%1.4d '];
end;
currStr = [currStr,'\n'];
fid=fopen(allBeta_file,'a');
fprintf(fid,'%d ',currIter);
fprintf(fid,currStr,betaOpt);
fclose(fid);

figure;
stem(betaOpt);
hold on;
stem(beta0,'r');
title(['betaOpt iter ',int2str(currIter),'vs previous iter (red)']);
%figName = ['Figs/betaOpt_',int2str(currIter),'_vsPrevious.fig'];
%saveas(h,figName);
