nbRdmPts = 100;
ptDim=51;
rdmPts = zeros(ptDim,nbRdmPts);
for k=1:nbRdmPts
    rdmPts(:,k) = get_rdmGT_uniformly('FminData/subNwk11c.mat');
end;

imagesc([allPts;rdmPts']);
%save('rdmPts100_unif.mat','rdmPts');
