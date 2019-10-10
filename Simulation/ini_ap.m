function [P_UL_ini,P_DL_ini,A_ini] = ini_ap(rho, fdcomm, radar)
% ini_ap performs the initialization for the alternating projection
% d_DL: number of downlink streams for each DL UE
% d_UL: number of uplink streams for each UL UE
%---- Intialize based on the channel matrices
I = fdcomm.DL_num;
J = fdcomm.UL_num;
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
    P_UL_i = P_UL(ii); %% Power level of the ith UL UE  
    HiB = H_UL{ii}; %% UL channel matrix with dimension Mc*Ni
    Ni     = size(HiB,2);
    d_iB = d_UL(ii);
    HiB_H = HiB';%HiB_H 
    P_iB_ini = HiB_H(:,1:d_iB); % Right singular matrix Ni*di
    [U_UL_i,~,V_UL_i] = svd(P_iB_ini);
    %s = diag(S);
    %[~,Idx] = sort(s,'descend');
%     U_sort = U(:,Idx);
%     V_sort = V(:,Idx);
    s_sort_nonzero = sqrt(P_UL_i/d_iB)*ones(d_iB,1);
    S_sort = [diag(s_sort_nonzero);zeros(Ni-d_iB,d_iB)];
    %S_sort = diag(s_sort_nonzero);
    P_iB_ini = U_UL_i*S_sort*V_UL_i';
    for k = 1 : K*N
        P_UL_ini{ii,k} = P_iB_ini;
    end
end
%% DL Precoding matrices
P_DL_ini = cell(J,K*N);
for jj = 1:J
    HBj = H_DL{jj}; %% UL channel matrix with dimension Ni*Mc
    Mc = size(HBj,2);
    d_Bj = d_DL(jj);
    HBj_H = HBj';
    P_Bj_ini = HBj_H(:,1:d_Bj);
    [U_DL_j,~,V_DL_j] = svd(P_Bj_ini);
    %s_DL = diag(S_DL);
    %[~,Idx] = sort(s_DL,'descend');
    %U_sort_j = U_DL(:,Idx);
    %V_sort_j = V_DL(:,Idx);
    s_sort_nonzero_j = sqrt(P_DL/J/d_Bj)*ones(d_Bj,1);
    S_sort_j = [diag(s_sort_nonzero_j);zeros(Mc-d_Bj,d_Bj)];
    P_Bj_ini = U_DL_j*S_sort_j*V_DL_j';
    for k = 1 : K*N
        P_DL_ini{jj,k} = P_Bj_ini;
    end
end
% Radar Code Matrix Intialization with PAR constraint
A_ini = zeros(K,Mr);
for mr = 1:Mr
    P_mr = P_r(mr);
    a_mr = P_mr/K;
    rho_mr = radar.rho(mr);
    a_mr_ini = Nearest_PAR(a_mr,rho_mr,P_mr);
    A_ini(:,mr) = a_mr_ini;
end
% Radar Code Matrix Intialization without PAR constraint
A_ini = zeros(K,Mr);
for mr = 1:Mr
    P_mr = P_r(mr);
    a_mr = P_mr/K;
    rho_mr = rho(mr);
    a_mr_ini = Nearest_PAR(a_mr,rho_mr,P_mr);
    A_ini(:,mr) = a_mr_ini;
end
end

