function [stress,R_iso,estrain,pstrain,peeq,num_iteration] = monotonic_hardening(mat_paras_mon,strain)
% This function is used to characterize uniaxial monotonic hardening behaviour

%% set tolerance
tol = 1.0e-8;
max_iter = 20;

%% monotonic plastic parameters
E_mod     = mat_paras_mon(1); % initial elastic modulus
sigma_y   = mat_paras_mon(2); % initial yield stress
strain_sh = mat_paras_mon(3); % length of yield plateau
sigma_sat = mat_paras_mon(4); % nonlinear hardening
m_N       = mat_paras_mon(5); % nonlinear hardening modulus
m_L       = mat_paras_mon(6); % linear hardening modulus

%% initialize state variables
% output variables
num_strain     = numel(strain);
stress         = zeros(num_strain,1);
R_iso          = ones(num_strain,1)*sigma_y;
estrain        = zeros(num_strain,1);
pstrain        = zeros(num_strain,1);
peeq           = zeros(num_strain,1);
num_iteration  = zeros(num_strain,1);

% temporary variables
flow_direction = zeros(num_strain,1);

%% calculate stress according to the strain history
for i = 2:numel(strain)
    
    % read the state variable from last step
    sigma_n      = stress(i-1);
    R_iso_n      = R_iso(i-1);
    estrain_n    = estrain(i-1);
    pstrain_n    = pstrain(i-1);
    peeq_n       = peeq(i-1);  
    np_n         = flow_direction(i-1); 
    
    % calculate strain increment
    d_strain_n = strain(i) - strain(i-1);
    
    % elastic prediction
    sigma_n1_trial = sigma_n + E_mod*d_strain_n;
    f_n1_trial     = abs(sigma_n1_trial) - R_iso_n;
    
    if f_n1_trial <= tol %elastic state
        
        % update output variables
        stress(i)        = sigma_n1_trial;
        R_iso(i)         = R_iso_n;
        estrain(i)       = estrain_n;
        pstrain(i)       = pstrain_n; 
        peeq(i)          = peeq_n;
        num_iteration(i) = 0;
        
        % update temporary variables
        flow_direction(i) = np_n;     
        
    else  % plastic correction
                  
        delta_gamma = 0.0;
        f_n1        = f_n1_trial;
        np_n1       = sign(sigma_n1_trial);
        num_iter    = 1;           
            
        while abs(f_n1) > tol && num_iter <= max_iter

            % calculate peeq
            peeq_n1 = peeq_n + delta_gamma;
            [R_iso_n1,dR_iso_n1] = mono_R(mat_paras_mon,peeq_n1);

            f_n1 = abs(sigma_n1_trial) - E_mod*delta_gamma - R_iso_n1;

            if abs(f_n1) > tol

                %iteration
                num_iter = num_iter + 1;
                
                if num_iter > max_iter
%                     'too many iterations'
                    break
                else
                    df_n1 = -E_mod - dR_iso_n1;
                    delta_gamma = delta_gamma - f_n1/df_n1;
                end

            end

        end

            % calculate state variables 
            sigma_n1 = sigma_n1_trial - E_mod*delta_gamma*np_n1;
            estrain_n1 = estrain_n + d_strain_n - delta_gamma*np_n1;
            pstrain_n1 = pstrain_n + delta_gamma*np_n1;
            
            % update output variables
            stress(i)         = sigma_n1;
            R_iso(i)          = R_iso_n1;
            estrain(i)        = estrain_n1;
            pstrain(i)        = pstrain_n1;
            peeq(i)           = peeq_n1;
            num_iteration(i)  = num_iter;

            % update temporary variables
            flow_direction(i)    = np_n1;
        
        end
        
end

end