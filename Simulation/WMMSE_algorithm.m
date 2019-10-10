function [fdcomm, radar, radar_comm] =...
    WMMSE_algorithm(fdcomm, radar, radar_comm, cov)
%-------------------------------------------------------------        
% WMMSE_algorithm returns the optimal precoding matrices ...
%   and radar transmit code matrix. 
% U_UL,U_DL,Ur,W_UL,W_DL,H_UL,H_DL,HrB,H_r_DL,H_UL_DL all cell arrays
% W_DL,Wr
K = radar.codelength;
I = fdcomm.UL_num;
J = fdcomm.DL_num;
iter_2_max = radar_comm.iter2_max;
% initialize precoding matrices and code matrix
iter_2 = 1;
% Xi_comm = zeros(K,1);
% for k = 1:K
%    Xi_comm(k) = fdcomm.Xi_UL(k) + fdcomm.Xi_DL(k);
% end
% Xi_mse = sum(Xi_comm)+ sum(radar.Xi_r);
% radar_temp = radar;
% fdcomm_temp = fdcomm;
while iter_2<= iter_2_max
    %% Alternating minimization    
    %%% Sub gradient method to find lambdas and mus
    % P_UL = fdcomm.ULpower;
    % P_DL = fdcomm.DLpower; 
    for k = 1: K
        %% UL precoder
        for ii = 1 : I
            [fdcomm,radar,cov] = UL_subgradient(fdcomm, radar_comm, radar,cov, ii,k);
        end
        %% DL precoder
        for jj = 1 : J
            [fdcomm,radar,cov] = DL_subgradient(fdcomm, radar_comm, radar,cov, ii,k);
        end
        %% radar code 
        radar = radar_code(fdcomm, radar, radar_comm, cov, k);
    end
    %%% Uodate the WMMSE receivers and weight matrices
    fdcomm = Comm_MMSE(fdcomm, cov);
    radar = radar_MMSE(radar, cov);
%     Xi_temp = 0;
%     for k = 1:K
%         fdcomm_temp = Xi_comm_k(fdcomm_temp, k);
%         Xi_temp = fdcomm_temp.Xi_UL(k)+ fdcomm_temp.Xi_DL(k) + Xi_temp;
%     end
%     radar_temp = Xi_radar(radar_temp);
%     Xi_temp = Xi_temp + sum(radar_temp.Xi_r);
    %% calculate the Xi
    iter_2 = iter_2 + 1;
end

