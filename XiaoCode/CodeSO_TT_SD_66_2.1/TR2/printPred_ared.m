function printPred_ared(currIter,mk_xk,mk_xTrial,f_sim_xk,f_sim_xTrial,add2file)

if currIter==1
    fid=fopen(add2file,'a');
    fprintf(fid,'iter mk_xk mk_xTrial f_sim_xk f_sim_xTrial\n');
    fclose(fid);
end;

fid=fopen(add2file,'a');
fprintf(fid,'%d %5.4f %5.4f %5.4f %5.4f\n',currIter,mk_xk,mk_xTrial,f_sim_xk,f_sim_xTrial);
fclose(fid);

%{
% check that identifiers are unique
fid=fopen(add2file,'r');
dataObs = textscan(fid,'%d %f');
fclose(fid);
if (length(unique(dataObs{1}))~=length(dataObs{1}))
    error('there are duplicate point IDs');
end;
%}