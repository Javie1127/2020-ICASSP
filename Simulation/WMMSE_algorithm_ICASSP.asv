function [fdcomm, radar, radar_comm,Xi_mse] =...
    WMMSE_algorithm_ICASSP(fdcomm, radar, radar_comm)
%-------------------------------------------------------------        
% WMMSE_algorithm updates the precoding matrices ...
%   and radar transmit code matrix.
%% Initialization
[fdcomm, radar] = ini_icassp(fdcomm, radar);
cov = covmat(fdcomm,radar,radar_comm);
fdcomm = Comm_MMSE(fdcomm,radar,cov);
radar = radar_MMSE(radar,cov);
fdcomm = Xi_comm_k(fdcomm,k);
radar = Xi_radar(radar);
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
Xi_mse = zeros(iter_2_max,1);
while iter_2<= iter_2_max
    disp(iter_2);
    %% Alternating minimization    
    %%% Sub gradient method to find lambdas and mus
    % P_UL = fdcomm.ULpower;
    % P_DL = fdcomm.DLpower; 
    for k = 1: K
        %% UL precoder
        for ii = 1 : I
            [fdcomm,radar,cov] = UL_subgradient(fdcomm, radar_comm, radar, cov, ii,k);
        end
        %% DL precoder
        for jj = 1 : J
            [fdcomm,radar,cov] = DL_subgradient(fdcomm, radar_comm, radar, cov, jj ,k);
        end
        %% radar code 
        [fdcomm, radar, cov] = radar_subgradient(fdcomm, radar_comm, radar, cov,k);
    end
    %%% Update the WMMSE receivers and weight matrices
    cov = covmat(fdcomm, radar, radar_comm);
    fdcomm = Comm_MMSE(fdcomm, radar, cov);
    radar = radar_MMSE(radar, cov);
    for k = 1:K
        fdcomm = Xi_comm_k(fdcomm, k);
        Xi_mse(iter_2) = fdcomm.Xi_UL(k)+ fdcomm.Xi_DL(k) + Xi_mse(iter_2);
    end
    radar = Xi_radar(radar);
    Xi_mse(iter_2) = Xi_mse(iter_2) + sum(radar.Xi_r);
    %% calculate the Xi
    iter_2 = iter_2 + 1;
end

