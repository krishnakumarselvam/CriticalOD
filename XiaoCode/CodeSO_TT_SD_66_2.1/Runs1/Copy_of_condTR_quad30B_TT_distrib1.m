% differs from condTR_quad30B_TT_distrib1 :
% - does not resort to queueing model

dbstop if error;


%%%%%
% 0. initializations
%%%%%

disp('for now we use xk values that are not rounded (denoted xStar), we round them to feed them to aimsun (gt_opts), this rounding might not be necessary');
disp('we use xStar values to estimate phi');

initIter = 4;
endIter = 5;
initPtId = 21;
directoryID = '1';

initPtS = int2str(initPtId);
% folderName = 'N:/RunsTR/Quad30/Runs',directoryID,'/'];
folderName = 'C:/MIT/Lab/KanchanaNanduri/CodeSO_TT_FC/Runs1/';
%folderName = '';

diary([folderName,'outMatlab_quad30B_TT_distrib.txt']);
diary on
dd=date;
cc=clock;
disp([dd,'  ',int2str(cc(4)),':',int2str(cc(5)),':',num2str(cc(6))]);

nbVblePhPerNode = [7,9,3,6,5,6,2,8,5];
csPh = cumsum(nbVblePhPerNode);
phases2Remove = [2,csPh(1)+9,csPh(2)+1,csPh(3)+3,csPh(4)+5,csPh(5)+3,csPh(6)+1,csPh(7)+8,csPh(8)+4];

for dbnIter = initIter:endIter
    
    iterDbnS = int2str(dbnIter);
    disp(['dbnIter : ',iterDbnS]);
    
    % Initialize some trust region parameters.
    %% initial values taken from MATLAB's trustnleqn.m script
    delta0 = 10;
    delta = delta0;
    deltaMax = 1e10;
    
    eta1 = .001;
    gammaInc = 1.2;
    gammaRed = .9;
    gammaMin = 10^-2;
    
    tolFmin = 10^-3; % fmincon tolerance for constraints and for rel. change in obj fctn value
    tolConstr = 10^-2; % fmincon tolerance for constraints and for rel. change in obj fctn value
    
    betaTol = 10^-1; % threshold for relative change in the norm of metamodel parameters (beta,alpha)
    mixturePtId = 4001; % Id of diversification points
    
    iter = 0;
    succIter = 0; % last successfull iteration
    stepAccept = 1;
    %iterMax = 50;
    sampleMax = 150;
    iterSample = 0;
    nbConsecFails = 0; % nb of consecutive trial points rejected
    nbFailsThresh = 10; % nb of falis threshold (reduces TR radius)
    tradeoff = 1.43;%value for SD
    runAimsun = 1; % so that we can progress on the script regardless of the 1001 issues with aimsun (license, console,..)
    
    fileNwk = ['../TR2/subNwk11c.mat'];
    qgParamsF = ['../TR2/qgParamsInit_cf4.mat'];
    nwkID = '11c';
    % xOut_k.mat represent x_out of the k^th iteration that will be used as the initial point (x_init) for the k+1^th iteration
    fileAdpatedPlans_txt = [folderName,'GTsUsed/oldVsNewPlans_quad30B_TT_',initPtS,'_',iterDbnS,'_0.txt'];
    file_gtOpt = [folderName,'GTsUsed/gtOpt_quad30B_TT_',initPtS,'_',iterDbnS,'_0.mat'];
    fileIn_py = [folderName,'PyInputs/new_gtDico_quad30B_TT_',initPtS,'_',iterDbnS,'_0.py'];
    %outFilePrint = 'FminData/outFmincon_0.txt'; % at iteration 0 its not used
    allObsFile = [folderName,'Obsvtns/allObs_quad30B_TT_',initPtS,'_',iterDbnS,'.txt'];
    allPtsFile = [folderName,'Obsvtns/allPoints_quad30B_TT_',initPtS,'_',iterDbnS,'.txt'];
    allBetaFile = [folderName,'Obsvtns/allBeta_quad30B_TT_',initPtS,'_',iterDbnS,'.txt'];
    rhoFile = [folderName,'Obsvtns/allPredAred_quad30B_TT_',initPtS,'_',iterDbnS,'.txt'];
    allRoundedGTs = [folderName,'Obsvtns/allRoundedGTs_quad30B_TT_',initPtS,'_',iterDbnS,'.txt'];
    allRelBetaF = [folderName,'Obsvtns/allRelBeta_quad30B_TT_',initPtS,'_',iterDbnS,'.txt'];
    tradeoff_py = [folderName,'PyInputs/tradeoffSD.py'];
    %addpath('C:/EPFL/QBD/QBDPiSqFullSyst/Runs/'); addpath('C:/EPFL/QBD/QBDPiSqFullSyst/');
    % addpath('N:/RunsTR/Quad30/QBD/'); addpath('N:/RunsTR/Quad30/TR2/'); addpath('.');
    addpath('../QBD/'); addpath('../TR2/'); addpath('.');
    printtradeoff(tradeoff_py,tradeoff);
    
    % used when get_fsim is called for a design point where x_out isnot yet available
    loadGtsFromQg_bool = 0;
    loadX0_inF = '';
    
    load('../TR2/rdmPts100_unif.mat','rdmPts');
    currGT_notFromQg = rdmPts(:,initPtId)';
    x_out = currGT_notFromQg;
    save(['FminData/xOut_quad30B_TT_',initPtS,'_',iterDbnS,'_',int2str(iter),'.mat'],'x_out');
    clear rdmPts;
    
    % calculate and print fSim(x_k)
    iterSample = iterSample+1;
    [x_k, f_sim_xk] = get_fSim1_TT(loadX0_inF,fileNwk,qgParamsF,fileAdpatedPlans_txt,fileIn_py,file_gtOpt,iter,allObsFile,allPtsFile,runAimsun,loadGtsFromQg_bool,currGT_notFromQg,allRoundedGTs,dbnIter,initPtId,folderName);
    
    if x_k ~=currGT_notFromQg error(''); end;
    
    ptDim = length(x_k);
    
    betaFile = estimatePoly_lsq_quad_dbn_augData8B(ptDim,allPtsFile,allObsFile,succIter,iter,allBetaFile,iterDbnS,'',initPtS,folderName,fileNwk,phases2Remove);
    load(betaFile,'betaOpt','aIndx','bIndx','cIndx');
    xkPoly = x_k;
    xkPoly(phases2Remove) = [];
    mk_xk = (betaOpt(aIndx)+xkPoly*betaOpt(bIndx)+(xkPoly.^2)*betaOpt(cIndx));
    clear betaOpt aIndx bIndx cIndx;
    
    disp(' ------- iterations beginning');
    
    % Beginning of main iteration loop.
    while iterSample <= sampleMax
        
        iter = iter + 1;
        iterS = int2str(iter);
        
        %%%%%
        % 1. Final criticality test : not used as termination criterion
        %%%%%
        
        %%%%%
        % 2. step calculation
        %%%%%
        % Compute step d
        saveX0_outF =  [folderName,'FminData/xOut_quad30B_TT_',initPtS,'_',iterDbnS,'_',iterS,'.mat'];
        outFilePrint =  [folderName,'FminData/outFmincon_quad30B_TT_',initPtS,'_',iterDbnS,'_',iterS,'.txt'];
        loadX0_inF =  [folderName,'FminData/xOut_quad30B_TT_',initPtS,'_',iterDbnS,'_',int2str(iter-1),'.mat'];
        
        [x_trial, mk_xTrial] = fmincon_totTimeGamma_someMuExog_cf6_Phi_quad30B_i(x_k,loadX0_inF,qgParamsF,fileNwk,delta,betaFile,outFilePrint,saveX0_outF,tolFmin,tolConstr,phases2Remove);
        
        %%%%%
        % 3. evaluation of trial point
        %%%%%
        
        % predicted model reduction given by d
        pred =  mk_xk - mk_xTrial; % calculate analytically (simplifications appear)
        
        % actual reduction
        % f_sim(x_trial): evaluate obj. fctn at trial point.
        fileAdpatedPlans_txt =  [folderName,'GTsUsed/oldVsNewPlans_quad30B_TT_',initPtS,'_',iterDbnS,'_',iterS,'.txt'];
        fileIn_py =  [folderName,'PyInputs/new_gtDico_quad30B_TT_',initPtS,'_',iterDbnS,'_',iterS,'.py'];
        file_gtOpt =  [folderName,'GTsUsed/gtOpt_quad30B_TT_',initPtS,'_',iterDbnS,'_',iterS,'.mat'];
        outFilePrint =  [folderName,'FminData/outFmincon_quad30B_TT_',initPtS,'_',iterDbnS,'_',iterS,'.txt'];
        loadX0_inF = saveX0_outF;
        
        iterSample = iterSample+1;
        [x_Trial_tmp, f_sim_xTrial] = get_fSim1_TT(loadX0_inF,fileNwk,qgParamsF,fileAdpatedPlans_txt,fileIn_py,file_gtOpt,iter,allObsFile,allPtsFile,runAimsun,loadGtsFromQg_bool,x_trial,allRoundedGTs,dbnIter,initPtId,folderName);
        
        if (x_trial~=x_Trial_tmp) error('inconsistent trial gts'); end;
        
        ared = f_sim_xk - f_sim_xTrial;
        
        rho = ared/pred;
        
        printPred_ared(iter,mk_xk,mk_xTrial,f_sim_xk,f_sim_xTrial,rhoFile);
        
        % define new iterate
        if (ared>0) % safeguard because fmincon can yield an optimal that has a larger obj fctn value than the initial point, because search direction is small and constraints are satisfied
            if ((rho > eta1))
                % accept step.
                stepAccept = 1;
                x_k = x_trial;
                f_sim_xk = f_sim_xTrial;
                %mk_xk = mk_xTrial; done later in step mdl update (which is done at each iteration since we include the new observation from fSim into the estimation of phi)
                succIter = iter;
                loadX0_inF = [folderName,'FminData/xOut_quad30B_TT_',initPtS,'_',iterDbnS,'_',iterS,'.mat'];
                nbConsecFails = 0;
                
            else stepAccept = 0;
                nbConsecFails = nbConsecFails+1;
            end;
        else stepAccept = 0;
            nbConsecFails = nbConsecFails+1;
        end;
        
        %%%%%
        % 4. Model improvement and update
        %%%%%
        
        % update model
        betaFile_prev = betaFile;
        betaFile = estimatePoly_lsq_quad_dbn_augData8B(ptDim,allPtsFile,allObsFile,succIter,iter,allBetaFile,iterDbnS,betaFile_prev,initPtS,folderName,fileNwk,phases2Remove);
        
        load(betaFile,'betaOpt','aIndx','bIndx','cIndx');
        xkPoly = x_k;
        xkPoly(phases2Remove) = [];
        mk_xk = (betaOpt(aIndx)+xkPoly*betaOpt(bIndx)+(xkPoly.^2)*betaOpt(cIndx));
        betaCurr = betaOpt;
        clear betaOpt aIndx bIndx cIndx;
        
        % compare succesive beta's
        % if necessary perform a diversification step
        load(betaFile_prev,'betaOpt');
        betaPrev = betaOpt;
        clear betaOpt;
        relBeta = norm(betaCurr-betaPrev)/norm(betaPrev);
        
        disp(['relBeta :',num2str(relBeta)]);
        printRelBeta(relBeta,allRelBetaF,iter);
        
        if (relBeta < betaTol) % sample another point.
            
            % retrieve info from exiting files, print into current outputfiles
            iterSample = iterSample+1;
            
            xOutF =  [folderName,'FminData/xOut_quad30B_TT_',initPtS,'_',iterDbnS,'_',int2str(mixturePtId),'.mat'];
            fileAdpatedPlans_txt =  [folderName,'GTsUsed/oldVsNewPlans_quad30B_TT_',initPtS,'_',iterDbnS,'_',int2str(mixturePtId),'.txt'];
            fileIn_py =  [folderName,'PyInputs/new_gtDico_quad30B_TT_',initPtS,'_',iterDbnS,'_',int2str(mixturePtId),'.py'];
            file_gtOpt =  [folderName,'GTsUsed/gtOpt_quad30B_TT_',initPtS,'_',iterDbnS,'_',int2str(mixturePtId),'.mat'];
            
            currGT_notFromQg = get_rdmGT_uniformly(fileNwk);
            
            % calculate and print fSim(x_k)
            get_fSim1_TT('',fileNwk,qgParamsF,fileAdpatedPlans_txt,fileIn_py,file_gtOpt,mixturePtId,allObsFile,allPtsFile,runAimsun,loadGtsFromQg_bool,currGT_notFromQg,allRoundedGTs,dbnIter,initPtId,folderName);
            
            % update mk
            betaFile_prev = betaFile;
            betaFile = estimatePoly_lsq_quad_dbn_augData8B(ptDim,allPtsFile,allObsFile,succIter,mixturePtId,allBetaFile,iterDbnS,betaFile_prev,initPtS,folderName,fileNwk,phases2Remove);
            
            mk_xkPrev = mk_xk; % so that it is stored in quad*.mat file
            load(betaFile,'betaOpt','aIndx','bIndx','cIndx');
            xkPoly = x_k;
            xkPoly(phases2Remove) = [];
            mk_xk = (betaOpt(aIndx)+xkPoly*betaOpt(bIndx)+(xkPoly.^2)*betaOpt(cIndx));
            clear betaOpt aIndx bIndx cIndx;
            
            mixturePtId = mixturePtId+1;
        end;
        
        %%%%%
        % 5. update radius
        %%%%%
        
        % Update trust region radius.
        if (rho > eta1)
            delta = min(gammaInc*delta,deltaMax);
        else
            if nbConsecFails == nbFailsThresh
                nbConsecFails = 0;
                delta = max(gammaRed*delta,gammaMin);
                disp('delta reduced');
            end;
        end;
        if delta == gammaMin
            error('launch safeguard method');
        end;
        
        save([folderName,'FminData/quad30B_TT_',initPtS,'_',iterDbnS,'_',iterS,'.mat']);
        
        if (mod(iter,3)==0) close all; end;
    end;
end;
diary off
