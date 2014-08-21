function [fQ,fQ_TT,SDALL,exitflag]  = get_fQ_rdmPt10(fileNwk,qgParamsF,xOutFile,currGTs,currPtId,iterDbnS,dirName,tradeoff)

load(fileNwk,'gamma_s','mu_s','serv_s','buffer_s','pij_s','signalizedQs_s','nonSignalizedQs_s','availableGreenTimes_s','phase2Node','cycleTimes_s','initialGreenTimePerPhase_s','q2phases','freeFlowCap','fixedGreenTime');


% set mu_s value according to the design point currGTs (which are normalized divided by their cycleTime)
mu_updated = mu_s;
for j =1:length(signalizedQs_s)
    q = signalizedQs_s(j);
    
    currPhases = q2phases{q};
    currNodeIndx = phase2Node(currPhases);
    currNodeIndx = unique(currNodeIndx);
    if (length(currNodeIndx)~=1) error('all phases should be associated to the same node'); end;
    % normalized vsn:
    % eqtn: mu := sum(green times of phases)*freeFlowCap + fixedGreenTime*freeFlowCap/cycleTime
    
    mu_updated(q) = sum(currGTs(currPhases))*freeFlowCap + fixedGreenTime(q)*freeFlowCap./cycleTimes_s(currNodeIndx);
end;

loadF='';

knownX0 = 0; % i.e. P_bs matrix can be loaded instead of recalculated
load_x0File = ''; %'x0_LIN_simLo.mat';

%%%%%%%%%%%%%%%% 1st run - initSet7
precision=1e-07;
myOptions{1}=precision;
prec = '7';

medAlgo{1}='';
myOptions{2}='trust-region-dogleg'; %'trust-region-reflective'; % 'on'; % medium algo == off
myOptions{3}=1000; % max iter
myOptions{4}='off'; % derivative check


iterR = 0;
iterLocal = 1;
exitflag=-1;
while ((iterLocal <= 2) & (exitflag<=0))
    
    switch iterLocal
        case 1
            scenarios{1} ='piUnif'; outPi='Unif';
            scenarios{2} =  'DeducedEqtn'; outMu='DedEqtn';
        otherwise
            scenarios{1} ='piRdm'; outPi='Rdm';
            rand('state', sum(100*clock));
            scenarios{2} =  'DeducedEqtn'; outMu='DedEqtn';
    end;
    %scenarios{2} =  'Mu'; outMu='Mu';
    
    outFile_lamFC = [dirName,'FminData/algolamFC_fQInit_fsolve_TT_11c_cf5_rdmPt_pi',outPi,'_mut',outMu,'_e0',prec,'_',iterDbnS,'_',int2str(currPtId),'.tab'];
    
    disp('starting lambda_FC calculation');
    lambda_flowCvtn = lambdaFlowCvtn_linear(gamma_s, pij_s,outFile_lamFC);
    disp('starting jacob calculation');
    
    outFile_pi_a = [dirName,'FminData/algo_fQInit_fsolve_11c_TT_cf5_rdmPt_pi',outPi,'_mut',outMu,'Large',myOptions{2},'_e0',prec,'_',iterDbnS,'_',int2str(currPtId),'.tab'];
    matFile_a = [dirName,'FminData/data_fQInit_fsolve_11c_TT_cf5_rdmPt_pi',outPi,'_mut',outMu,'Large',myOptions{2},'_e0',prec,'_',iterDbnS,'_',int2str(currPtId)];
    
    tic
    [x_out,fval,exitflag,output,jacobian,extra_vbles_out] = nwk_fsolve_muHat_cf5(gamma_s, lambda_flowCvtn, mu_updated, serv_s, buffer_s, pij_s,outFile_pi_a,scenarios,myOptions,knownX0,load_x0File,iterR);
    tTime = toc;
    disp(['ended iterLocal: ',int2str(iterLocal)]);
    iterLocal = iterLocal +1;
end;

save(matFile_a);
disp('initial run: adapting xOut from fsolve to fmincon format');
x_out = [x_out,mu_updated(signalizedQs_s),currGTs];
save(xOutFile,'x_out');

if exitflag>0
    
    % calculate corresponding avgTT: fQ
    load(qgParamsF,'last_index','nb','nbVbles');
    load(fileNwk,'cycleTimes_s','condFact','phase2Node','gamma_s','buffer_s');
    nbPhases = length(phase2Node);
    gamma = gamma_s;
    buffer = buffer_s;
    clear gamma_s buffer_s;
    xFull = x_out;
    clear x_out;
    b=max(buffer)+1;
    otmp=0:1:max(buffer);
    om=repmat(otmp,nb,1);
    
    Inside=zeros(nb,b);
    sumInside=zeros(nb,b);
    for p=1:nb;
        for o = 0:buffer(p);
            Inside(p,o+1)=(2*o+2);
        end;
    end;
    for p=1:nb;
        for o = 0:buffer(p);
            sumInside(p,o+1)=sum(Inside(p,1:o+1));
        end;
    end;
    
    rho = xFull(nb*2+(1:nb))./xFull(nb+(1:nb));
    allNi = rho.*(1./(1-rho)-( ((buffer+2).*(rho.^(buffer+1)))./(1-rho.^(buffer+2)) ));
    %allNi_z2 = sum(allNi);
    %hasGammaInd_z2 = find(gamma>0);
    
    
    allG_zsd = xFull(nb*2+(1:nb)).*(1-xFull(nb*3 + (1:nb)));
    %allG_z2 = sum(gamma(hasGammaInd_z2).*(1-xFull(nb*3 + hasGammaInd_z2)));
    
    %fQ_TT = 3600.*(allNi_z2/(condFact*allG_z2));
    EXSD=allNi./(allG_zsd*(condFact/3600));
    fQ_TT = sum(EXSD);
    E2X = EXSD.^2;
    sumOut=zeros(nb,b);
    rhotmp=repmat(rho',1,b);
    
    sumOut=rhotmp.^om.*sumInside;
    
    EX2=((1./xFull(nb+(1:nb)))*1/(condFact/3600))'.^2.*(1-rho')./(1-rho'.^(buffer'+1)).*sum(sumOut(:,:),2);
    
    VAR=EX2' - E2X;
    VARALL=sum(VAR);
    SDALL=VARALL^0.5;
    
    fQ = fQ_TT+tradeoff*SDALL;
    
    clear x_out fval output jacobian pi_out extra_vbles_out gamma_s mu_s serv_s buffer_s pij_s signalizedQs_s nonSignalizedQs_s availableGreenTimes_s phase2Node cycleTimes_s initialGreenTimePerPhase_s q2phases lambda_flowCvtn;
else
    fQ=-1;
end;
