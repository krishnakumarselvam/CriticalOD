% This is a temporary script only. It reads from a temporary mat file the
% gamma, mu and pij, and solves the 6 variable queuing equation.

clc
clear
load Large.mat mu x_out nb_Pij nb nb_paths nb_OD Pij_rowindex Pij_colindex buffer serv gamma_index 

% Getting all the indices
% ######################################################################################################################
%%
for i=1:1
Rho_hat_index = 1:nb;
Rho_eff_index = nb+1:2*nb;
Lambda_eff_index = 2*nb+1:3*nb;
Pnk_index=3*nb+1:4*nb;
Pi_full_index = 4*nb+1:5*nb;
ENi_index=5*nb+1:6*nb;
Ti_index=6*nb+1:7*nb;
Ci_index=7*nb+1:7*nb+nb_paths;
OD_sum_index=7*nb+nb_paths+1:7*nb+nb_paths+nb_OD;
Lij_index=7*nb+nb_paths+nb_OD+1:7*nb+nb_paths+nb_OD+nb_paths;
F_index=7*nb+nb_paths+nb_OD+nb_paths+1:7*nb+nb_paths+nb_OD+nb_paths+nb_paths;
gamma_index=7*nb+nb_paths+nb_OD+nb_paths+nb_paths+1:7*nb+nb_paths+nb_OD+nb_paths+nb_paths+nb;
Q_Pij_index=7*nb+nb_paths+nb_OD+nb_paths+nb_paths+nb+1:7*nb+nb_paths+nb_OD+nb_paths+nb_paths+nb+nb_Pij;
true_rho_index=7*nb+nb_paths+nb_OD+nb_paths+nb_paths+nb+nb_Pij+1:7*nb+nb_paths+nb_OD+nb_paths+nb_paths+nb+nb_Pij+nb;
end
% ######################################################################################################################
%%
%Load relevant data
 
gamma = x_out(gamma_index);
mu = mu; 
Pij = sparse(Pij_rowindex,Pij_colindex,x_out(Q_Pij_index),nb,nb);

[solution Exp_length Delay_data] = QueuingFsolve(gamma, Pij, buffer,serv,mu);