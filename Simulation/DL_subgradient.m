function [fdcomm, radar, cov] = DL_subgradient(fdcomm, radar_comm, radar, cov, jj, k)
% sub-graident method to find 
% update only the P_dj_k
I = fdcomm.UL_num;
J = fdcomm.DL_num;
t_DL_max = fdcomm.t_DL_max;
P_B = fdcomm.DLpower;
t_DL = 1;
R_DL = fdcomm.R_DL; % UL minimum rate
fdcomm_temp = fdcomm; 
% radar_temp = radar;
Xi_k = fdcomm.Xi_UL(k)+fdcomm.Xi_DL(k) + sum(radar.Xi_r);
lambda_k_d_t = fdcomm.lambda_DL(k);
mu_jd_k_t = fdcomm.mu_DL(jj,k);
d_DL_j = fdcomm.DLstream_num(jj);
HBj = fdcomm.DLchannels{jj};
R_in_dj = cov.in_DL{jj,k};
while t_DL <= t_DL_max
    beta_d_k_t = 1/t_DL;
    epsilon_j_k_t = 1/t_DL;
    % update lambda_k_t
    P_jd_k_sum = 0;
    for jj = 1:J
        P_jd_k = fdcomm_temp.DLprecoders{jj,k};
        P_jd_k_sum = trace(P_jd_k*P_jd_k') + P_jd_k_sum;
    end
    lambda_k_d_temp = lambda_k_d_t + beta_d_k_t*(P_jd_k_sum-P_B);
    lambda_k_d_temp = max(0,lambda_k_d_temp);
    fdcomm_temp.lambda_DL(k) = lambda_k_d_temp;
    % update mu_j_k_t
%     E_jd_k_t = fdcomm.DL_MMSE{jj,k};
%     R_jd_k_t = log2(det((E_jd_k_t)^(-1)));
    R_jd_k_t = log2(det(eye(d_DL_j)+P_jd_k'*HBj'/R_in_dj*HBj*P_jd_k));
    mu_jd_k_temp = mu_jd_k_t + epsilon_j_k_t*(R_DL-R_jd_k_t);
    mu_jd_k_temp = max(mu_jd_k_temp,0);
    fdcomm_temp.mu_DL(jj,k) = mu_jd_k_temp;
    % update P_jd_k
    fdcomm_temp = DL_precoders(k, jj, fdcomm_temp, radar_comm, radar, cov);
%     % update cov_temp
%     cov_temp = covmat(fdcomm_temp, radar, radar_comm);
%     % update MMSE
     fdcomm_temp = Comm_MMSE(fdcomm_temp, radar, cov);
%     radar_temp = radar_MMSE(radar_temp, cov_temp);
    % update Xi_MSE
    fdcomm_temp = Xi_comm_k(fdcomm_temp, k);
%     radar_temp = Xi_radar(radar_temp);
    Xi_temp = fdcomm_temp.Xi_UL(k)+ fdcomm_temp.Xi_DL(k) + sum(radar.Xi_r);
    if Xi_temp < Xi_k
        fdcomm = fdcomm_temp;
        Xi_k = Xi_temp;
    end
    t_DL = t_DL+1;
end
end

