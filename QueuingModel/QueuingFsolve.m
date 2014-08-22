function [solution Exp_length Delay_data] = QueuingFsolve(gamma_s, pij_s, buffer_s,serv_s,mu_s)

nb=length(buffer_s);
mu_updated = mu_s;

knownX0 = 0; % i.e. P_bs matrix can be loaded instead of recalculated
load_x0File = ''; 

precision=1e-07;
myOptions{1}=precision;
myOptions{2}='trust-region-dogleg';
myOptions{3}=10000; % max iter
myOptions{4}='off'; % derivative check
iterR = 0;

scenarios{1} ='piUnif'; 
scenarios{2} =  'DeducedEqtn'; 


lambda_flowCvtn = Q_lamFC(gamma_s, pij_s);

[x_out,fval,exitflag,output,jacobian,extra_vbles_out] = SolveEquations(gamma_s, lambda_flowCvtn, mu_updated, serv_s, buffer_s, pij_s,scenarios,myOptions,knownX0,load_x0File,iterR);
disp('Queuing Model Solved');
rho_eff_index=nb*5;
lambda_eff_index=2*nb;
Pnk_index=3*nb;
solution=x_out';
sol_lambda_effective=solution(lambda_eff_index+1:lambda_eff_index+nb,1);
sol_Pni_ki=solution(Pnk_index+1:Pnk_index+nb,1);
sol_rho=solution(rho_eff_index+1:rho_eff_index+nb,1);
Exp_length=zeros(nb,1);
Delay=zeros(nb,1);
Total_Delay_Nr=0;
Total_Delay_Dr=0;

for i=1:nb
    rho=sol_rho(i,1);
    sol_k=buffer_s(i);
    Exp_length(i,1)=rho*(1/(1-rho) - ((1 + sol_k)*rho^sol_k)/(1-rho^(sol_k+1)) );
    Delay(i,1) = Exp_length(i,1)/ (sol_lambda_effective(i,1))  +  (sol_k- Exp_length(i,1))*4/60000;
    Total_Delay_Nr=Total_Delay_Nr+Exp_length(i,1);
    Total_Delay_Dr=Total_Delay_Dr+gamma(i)*(1-sol_Pni_ki(i,1));
end

Delay_data(1,1:nb)=Delay;
Delay_data(isnan(Delay_data))=0;
solution=solution';
Exp_length=Exp_length';

end