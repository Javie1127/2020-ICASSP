function [fdcomm, radar] = ini_icassp(fdcomm, radar)
% ini_ap performs the initialization for the alternating projection
% d_DL: number of downlink streams for each DL UE
% d_UL: number of uplink streams for each UL UE
%---- Intialize based on the channel matrices
I = fdcomm.UL_num;
J = fdcomm.DL_num;
P_UL = fdcomm.ULpower;
P_DL = fdcomm.DLpower;
Mr = radar.TX;
K = radar.codelength;
N = radar.PRI_num;
H_UL = fdcomm.ULchannels;
H_DL = fdcomm.DLchannels;
d_DL = fdcomm.DLstream_num;
d_UL = fdcomm.ULstream_num;
P_r = radar.power;
%% Strategy 1 Right singular value of Channel Matrices
% UL Precoding Matrices  
P_UL_ini = cell(I,K*N);
for ii = 1:I
    P_UL_i  = P_UL(ii); %% Power level of the ith UL UE  
    HiB     = H_UL{ii}; %% UL channel matrix with dimension Mc*Ni
    d_ui    = d_UL(ii);
    [~,~,V] = svd(HiB);
    P_iB_ini = V(:,1:d_ui);
    [U_UL_i,s,V_UL_i] = svd(P_iB_ini);
    s(1:d_ui,:) = diag(sqrt(P_UL_i/d_ui)*ones(d_ui,1));
    %S_sort = diag(s_sort_nonzero);
    P_iB_ini = U_UL_i*s*V_UL_i';
    for k = 1 : K
        P_UL_ini{ii,k} = P_iB_ini;
    end
end
% for ii = 1:I
%     P_UL_i  = P_UL(ii); %% Power level of the ith UL UE  
%     d_ui    = d_UL(ii);
%     P_iB_ini = 
%     for k = 1 : K*N
%         P_UL_ini{ii,k} = P_iB_ini;
%     end
% end
%% DL Precoding matrices
P_DL_ini = cell(J,K*N);
for jj = 1:J
    HBj = H_DL{jj}; %% UL channel matrix with dimension Ni*Mc
    %Mc = size(HBj,2);
    d_Bj = d_DL(jj);
    [~,~,V] = svd(HBj);
    P_dj_ini = V(:,1:d_Bj);
    [U_DL_j,s,V_DL_j] = svd(P_dj_ini);
    s(1:d_Bj,:) = diag(sqrt(P_DL/J/d_Bj)*ones(d_Bj,1));
    %HBj_H = HBj';
    %P_Bj_ini = HBj_H(:,1:d_Bj);
    %[U_DL_j,~,V_DL_j] = svd(P_Bj_ini);
    %s_DL = diag(S_DL);
    %[~,Idx] = sort(s_DL,'descend');
    %U_sort_j = U_DL(:,Idx);
    %V_sort_j = V_DL(:,Idx);
    %s_sort_nonzero_j = sqrt(P_DL/J/d_Bj)*ones(d_Bj,1);
    %S_sort_j = [diag(s_sort_nonzero_j);zeros(Mc-d_Bj,d_Bj)];
%     P_Bj_ini = U_DL_j*S_sort_j*V_DL_j';
    P_Bj_ini = U_DL_j*s*V_DL_j';
    for k = 1 : K
        P_DL_ini{jj,k} = P_Bj_ini;
    end
end

% Initialize as the uncoded waveform

A_ini = zeros(K,Mr);
for k = 1:K
    P_k = P_r(k);
    A_ini(:,k) = sqrt(P_k/Mr);
end
radar.codematrix                = A_ini;
fdcomm.ULprecoders              = P_UL_ini;
fdcomm.DLprecoders              = P_DL_ini;
fdcomm.R_UL = 5; 
fdcomm.R_DL = 10;
%% Subgradient method 
fdcomm.lambda_UL = ones(I,K);
fdcomm.lambda_DL = ones(K,1);
fdcomm.mu_UL = ones(I,K);
fdcomm.mu_DL = ones(J,K);
fdcomm.t_UL_max = 8; % max number of iterations to execute the UL subgradient method
fdcomm.t_DL_max = 8; % max number of iterations to execute the DL suggradient method
radar.t_r_max = 8;
radar.lambda = ones(K,1);
%% WMMSE 
fdcomm.Xi_UL = zeros(K,1);
fdcomm.Xi_DL = zeros(K,1);
end

