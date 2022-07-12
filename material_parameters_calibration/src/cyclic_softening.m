function [stress,alpha,Y_iso,estrain,pstrain,peeq,psi,phi,E_modulus,num_iteration] = ...
         cyclic_softening(mat_paras,num_alpha,num_Y_iso,strain)

% This function is used to characterize cyclic softening behaviour within the
% framework of metallic cyclic plasticity
% It can be used to calculate the stress response of with resepect to a given 
% strain history.

%% set tolerance
tol = 1.0e-8;
max_iter = 20;

%% assign material parameters
% mat_paras is a column vector that stores material parameters.
% elastic parameters
E_mod0       = mat_paras(1); % initial elastic modulus

% monotonic plastic parameters
sigma_y      = mat_paras(2); % initial yield stress
strain_sh    = mat_paras(3); % length of yield plateau
sigma_sat    = mat_paras(4); % nonlinear hardening saturation
m_N          = mat_paras(5); % nonlinear hardening modulus
m_L          = mat_paras(6); % linear hardening modulus

% cyclic plastic parameters
m_phi        = mat_paras(7);
phi_sat      = mat_paras(8);
m_alpha_list = mat_paras(9:8+num_alpha);
omega_list   = mat_paras(9+num_alpha:8+2*num_alpha);
m_y_list     = mat_paras(9+2*num_alpha:8+2*num_alpha+num_Y_iso);
Q_list       = mat_paras(9+2*num_alpha+num_Y_iso:8+2*num_alpha+2*num_Y_iso);
E_mod_sat    = mat_paras(9+2*num_alpha+2*num_Y_iso);
xi_E         = mat_paras(10+2*num_alpha+2*num_Y_iso);

%% initialize state variables
% output variables
num_strain     = numel(strain);
stress         = zeros(num_strain,1);
alpha          = zeros(num_strain,1);
Y_iso          = ones(num_strain,1)*sigma_y;
estrain        = zeros(num_strain,1);
pstrain        = zeros(num_strain,1);
peeq           = zeros(num_strain,1);
psi            = zeros(num_strain,1);
phi            = zeros(num_strain,1);
alpha_list     = zeros(num_strain,num_alpha);
Y_iso_list     = zeros(num_strain,num_Y_iso);
E_modulus      = ones(num_strain,1)*E_mod0;
num_iteration  = zeros(num_strain,1);

% temporary variables
flow_direction = zeros(num_strain,1);
g_emodulus     = ones(num_strain,1);

% state variables for monotonic load mode
state_variables_n = zeros(1,num_alpha);

%% calculate stress according to the strain history
for i = 2:numel(strain)
    
    % read the state variable from last step
    sigma_n      = stress(i-1);
    alpha_n      = alpha(i-1);
    Y_iso_n      = Y_iso(i-1);
    estrain_n    = estrain(i-1);
    pstrain_n    = pstrain(i-1);
    peeq_n       = peeq(i-1);
    psi_n        = psi(i-1);   
    phi_n        = phi(i-1);
    alpha_list_n = alpha_list(i-1,:);
    Y_iso_list_n = Y_iso_list(i-1,:);  
    np_n         = flow_direction(i-1); 
    g_n          = g_emodulus(i-1);
    
    % calculate strain increment
    d_strain_n = strain(i) - strain(i-1);
    
    % elastic prediction
    sigma_n1_trial = sigma_n + g_n*E_mod0*d_strain_n;
    xi_trial       = sigma_n1_trial - alpha_n;
    f_n1_trial     = abs(xi_trial) - Y_iso_n;
    
    if f_n1_trial <= tol %elastic state
        
        % update output variables
        stress(i)         = sigma_n1_trial;
        alpha(i)          = alpha_n;
        Y_iso(i)          = Y_iso_n;
        estrain(i)        = estrain_n;
        pstrain(i)        = pstrain_n; 
        peeq(i)           = peeq_n;
        psi(i)            = psi_n;
        phi(i)            = phi_n;
        alpha_list(i,:)   = alpha_list_n; 
        Y_iso_list(i,:)   = Y_iso_list_n;
        E_modulus(i)      = g_n*E_mod0;
        num_iteration(i)  = 0;
        
        % update temporary variables
        flow_direction(i) = np_n;
        g_emodulus(i)     = g_n;
        
        
    else  % plastic correction
                  
        delta_gamma = 0.0;
        f_n1        = f_n1_trial;
        num_iter    = 1;           
            
        while abs(f_n1) > tol && num_iter <= max_iter

            % calculate peeq
            peeq_n1 = peeq_n + delta_gamma;
            
            % calculate g_emodulus, xi_bar and flow direction
            g_n1   = (g_n + xi_E*delta_gamma*E_mod_sat/E_mod0)/(1.0 + xi_E*delta_gamma);
            dg_n1  = xi_E*(E_mod_sat/E_mod0 - g_n)*delta_gamma/(1.0 + xi_E*delta_gamma)^2;
            
            xi_bar = g_n1/g_n*sigma_n1_trial - sum(alpha_list_n./(1.0 + m_alpha_list*delta_gamma));
            np_n1  = sign(xi_bar);
            
            % calculate Y_iso_n1 and dY_iso_n1
            Y_iso_list_n1  = (Y_iso_list_n + m_y_list.*Q_list*delta_gamma)./(1.0 + m_y_list*delta_gamma);
            dY_iso_list_n1 = m_y_list.*(Q_list - Y_iso_list_n)./(1.0 + m_y_list*delta_gamma).^2;
            Y_iso_n1       = sigma_y + sum(Y_iso_list_n1);
            dY_iso_n1      = sum(dY_iso_list_n1);
            
            % calculate psi_n1           
            if psi_n == 0 && np_n*np_n1 < 0.0
                psi_n1 = 1;
            else
                psi_n1 = psi_n;
            end
            
            % calculate phi_n1 and dphi_n1
            if psi_n1 == 0
                if delta_gamma == 0.0
                    phi_n1  = phi_n;
                    dphi_n1 = 0.0;
                else
                    [phi_n1,dphi_n1,state_variables_n1] = phi_calc(mat_paras,state_variables_n,Y_iso_n1,...
                                                                   dY_iso_n1,peeq_n1,delta_gamma,num_alpha);
                end
            else
                phi_n1  = (phi_n + m_phi*phi_sat*delta_gamma)/(1.0 + m_phi*delta_gamma);
                dphi_n1 = m_phi*(phi_sat - phi_n)/(1.0 + m_phi*delta_gamma)^2;
            end
            

            f_n1 = abs(xi_bar) - (g_n1*E_mod0 + sum(m_alpha_list.*omega_list*phi_n1./(1.0 + ...
                   m_alpha_list*delta_gamma)))*delta_gamma - Y_iso_n1;

            if abs(f_n1) > tol

                %iteration
                num_iter       = num_iter + 1;
                
                if num_iter > max_iter
%                     'too many iterations'
                    break
                else
                    dxi_bar_dgamma = np_n1*(sigma_n1_trial/g_n*dg_n1 + ...
                                     sum(m_alpha_list.*alpha_list_n./(1.0 + m_alpha_list*delta_gamma).^2));

                    df_n1          = dxi_bar_dgamma - g_n1*E_mod0 - E_mod0*delta_gamma*dg_n1 - ...
                                     sum(m_alpha_list.*omega_list*phi_n1./(1.0 + m_alpha_list*delta_gamma).^2) - ...
                                     sum(m_alpha_list.*omega_list*delta_gamma./(1.0 + m_alpha_list*delta_gamma))*dphi_n1 - dY_iso_n1;

                    delta_gamma = delta_gamma - f_n1/df_n1;
                end

            end

        end

            % calculate state variables 
            sigma_n1 = g_n1/g_n*sigma_n1_trial - g_n1*E_mod0*delta_gamma*np_n1;
            alpha_list_n1 = (alpha_list_n + m_alpha_list.*omega_list*phi_n1*delta_gamma*np_n1)./(1.0 + m_alpha_list*delta_gamma);
            alpha_n1 = sum(alpha_list_n1);
            estrain_n1 = estrain_n + d_strain_n - delta_gamma*np_n1;
            pstrain_n1 = pstrain_n + delta_gamma*np_n1;
            
            % update output variables
            stress(i)         = sigma_n1;
            alpha(i)          = alpha_n1;
            Y_iso(i)          = Y_iso_n1;
            estrain(i)        = estrain_n1;
            pstrain(i)        = pstrain_n1;
            peeq(i)           = peeq_n1;
            psi(i)            = psi_n1;
            phi(i)            = phi_n1;
            alpha_list(i,:)   = alpha_list_n1;
            Y_iso_list(i,:)   = Y_iso_list_n1;
            E_modulus(i)      = g_n1*E_mod0;
            num_iteration(i)  = num_iter;

            % update temporary variables
            flow_direction(i)    = np_n1;
            g_emodulus(i)        = g_n1;
            
            % update state variables
            state_variables_n = state_variables_n1;
        
        end
        
end
    
end
