function [phi_n1,dphi_n1,state_variables_n1] = phi_calc(mat_paras,state_variables_n,Y_iso_n1,dY_iso_n1,...
                                                        peeq_n1,delta_gamma,num_alpha)
% This function is used to calculate the phi_n1 and dphi_n1

%% assign material parameters
% mataParas is a column vector that stores material parameters.
% elastic parameters
E_mod0       = mat_paras(1); % initial elastic modulus

% monotonic plastic parameters
sigma_y      = mat_paras(2); % initial yield stress
strain_sh    = mat_paras(3); % length of yield plateau
sigma_sat    = mat_paras(4); % nonlinear hardening saturation
m_N          = mat_paras(5); % nonlinear hardening modulus
m_L          = mat_paras(6); % linear hardening modulus

% cyclic plastic parameters
m_alpha_list = mat_paras(9:8+num_alpha);
omega_list   = mat_paras(9+num_alpha:8+2*num_alpha);

%% state variables
alpha_list_n = state_variables_n(1:num_alpha);

%% calculate monotonic hardening
[R_n1,dR_n1] = mono_R(mat_paras,peeq_n1);

%% calculate phi_m and dphi_m_n1
G         = R_n1 - Y_iso_n1 - sum(alpha_list_n./(1.0 + m_alpha_list*delta_gamma));
F         = sum(m_alpha_list.*omega_list*delta_gamma./(1.0 + m_alpha_list*delta_gamma));
dG_dgamma = dR_n1 - dY_iso_n1 + sum(m_alpha_list.*alpha_list_n./(1.0 + m_alpha_list*delta_gamma).^2);
dF_dgamma = sum(m_alpha_list.*omega_list./(1.0 + m_alpha_list*delta_gamma).^2);

phi_n1    = G/F;
dphi_n1   = (dG_dgamma*F - dF_dgamma*G)/F^2;

%% update alpha_list and alpha
alpha_list_n1 = (alpha_list_n + m_alpha_list.*omega_list*phi_n1*delta_gamma)./(1.0 + m_alpha_list*delta_gamma);

%% update state variables
state_variables_n1 = alpha_list_n1;

end
