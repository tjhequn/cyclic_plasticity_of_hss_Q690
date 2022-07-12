function [R_mon,dR_mon] = mono_R(mat_paras_mon,peeq)
% This function is used to compute monotonic hardening function R and its
% derivate.

%% monotonic plastic parameters
E_mod     = mat_paras_mon(1); % initial elastic modulus
sigma_y   = mat_paras_mon(2); % initial yield stress
strain_sh = mat_paras_mon(3); % length of yield plateau
sigma_sat = mat_paras_mon(4); % nonlinear hardening
m_N       = mat_paras_mon(5); % nonlinear hardening modulus
m_L       = mat_paras_mon(6); % linear hardening modulus

%% calculate monotonic hardening function R
if peeq > strain_sh % strain hardening stage
    temp_peeq = peeq - strain_sh;
    R_mon     = sigma_y + sigma_sat*(1.0 - exp(-m_N*temp_peeq)) + m_L*temp_peeq;
    dR_mon    = sigma_sat*m_N*exp(-m_N*temp_peeq) + m_L;
    
    R_mon = sigma_y*(1.0+strain_sh)/(1.0-sigma_y/E_mod)+sigma_sat*(1.0-exp(-m_N*temp_peeq))+m_L*temp_peeq;
    dR_mon = sigma_sat*m_N*exp(-m_N*temp_peeq)+m_L;

else % yield plateau
    R_mon  = sigma_y;
    dR_mon = 0.0;
    
    R_mon = sigma_y*(1.0+peeq)/(1.0-sigma_y/E_mod);
    dR_mon = sigma_y/(1.0-sigma_y/E_mod);

end

end