function [fdcomm, radar, cov] = radar_subgradient(fdcomm, radar_comm, radar, cov,k)
% UL subgradient method 
%
t_r_max = radar.t_r_max;
t_r = 1;
K = radar.codelength;
fdcomm_temp = fdcomm; 
radar_temp = radar;
cov_temp = cov;
fdcomm = Xi_comm_k(fdcomm,k);
radar = Xi_radar(radar);
Xi_k = fdcomm.Xi_UL(k)+fdcomm.Xi_DL(k) + sum(radar.Xi_r);
lambda_r_k_t = radar.lambda(k);
P_rk = sum(radar.power/K);
while t_r <= t_r_max
    % fdcomm tracks the optimal results
    % fdcomm_temp tracks the instantaneous lambda, Piu,
    beta_i_k_t = 1/sqrt(t_r); % step number
    a_k_t = radar_temp.codematrix(k,:).';
    % update lambda
    lambda_r_k_new = lambda_r_k_t + beta_i_k_t*(abs(a_k_t'*a_k_t)-P_rk);
    lambda_r_k_t = max(lambda_r_k_new,0);
    radar.lambda(k) = lambda_r_k_t;
    % Update a_k with lambda_r_k_t
    radar_temp = radar_code(fdcomm_temp, radar_temp, radar_comm, cov_temp,k);
    % Update covariance matrices 
    cov_temp = covmat(fdcomm_temp,radar_temp,radar_comm); 
    % update Comm WMMSE receivers and weight matrices
    fdcomm_temp = Comm_MMSE(fdcomm_temp, radar_temp, cov_temp);
    % Update radar WMMSE receivers and weight matrices
    radar_temp = radar_MMSE(radar_temp, cov_temp);
    % update Xi_MSE
    fdcomm_temp = Xi_comm_k(fdcomm_temp, k);
    radar_temp = Xi_radar(radar_temp);
    Xi_temp = fdcomm_temp.Xi_UL(k)+fdcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
    if Xi_temp < Xi_k
        fdcomm = fdcomm_temp;
        radar = radar_temp;
        cov = cov_temp;
        Xi_k = Xi_temp; 
    end
    % calculate the new Xi_mse
    t_r = t_r+1;
end
end

