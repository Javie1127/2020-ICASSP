function [fdcomm] = UL_precoders(k, ii, fdcomm, radar,cov)

I = fdcomm.UL_num;
J = fdcomm.DL_num;
% napla_PiB_R_UL_k = cell(I,1);
% napla_PiB_R_DL_k = cell(J,1);
napla_PiB_R_UL_sum = 0;
napla_PiB_R_DL_sum = 0;
d_UL = fdcomm.ULstream_num; % number of data streams of the UL UE
d_DL = fdcomm.DLstream_num;
H_iB = fdcomm.ULchannels{ii,1};
U_iu_k = fdcomm.UL_WMMSE_RX{ii,k};
W_iu_k = fdcomm.UL_weights{ii,k};
tilde_P_iu_k = fdcomm.ULprecoders{ii,k}; % being updated
inverse_UL = 0;
alpha_UL = fdcomm.alpha_UL;
alpha_iu = alpha_UL(ii);
lambda_iu_k = fdcomm.lambda_UL(ii,k);
%% Derivatives of the UL/DL achivable rate
for q = 1:I
    R_in_uq_k = cov.in_UL{q,k};
    H_qB = fdcomm.ULchannels{q,1};
    P_qu_k = fdcomm.ULprecoders{q,k};
    mu_qu_k = fdcomm.mu_UL(q,k);
    if q == ii
        napla_PiB_R_UL_sum = 2*mu_qu_k*H_qB'/R_in_uq_k*H_qB*P_qu_k/...
           (eye(d_UL(q))+P_qu_k'*H_qB'/R_in_uq_k*H_qB*P_qu_k) + napla_PiB_R_UL_sum;
    else
        napla_PiB_R_UL_sum = -2*mu_qu_k*H_qB'/R_in_uq_k*H_qB*P_qu_k/...
           (eye(d_UL(ii))+P_qu_k'*H_qB'/R_in_uq_k*H_qB*P_qu_k)*...
           P_qu_k'*H_qB'/R_in_uq_k*H_iB*tilde_P_iu_k+napla_PiB_R_UL_sum;
    end
    alpha_q = alpha_UL(q);
    U_qu_k = fdcomm.UL_WMMSE_RX{q,k};
    W_qu_k = fdcomm.UL_weights{q,k};
    inverse_UL = alpha_q*H_iB'*U_qu_k'*W_qu_k*U_qu_k*H_iB+inverse_UL;
end
inverse_DL = 0;
alpha_DL = fdcomm.alpha_DL;
for jj = 1:J
    mu_gd_k = fdcomm.mu_DL(jj,k);
    H_i_j = fdcomm.ULDLchannels{ii,jj};
    H_Bj = fdcomm.DLchannels{jj};
    R_in_d_j_k = cov.in_DL{jj,k};
    P_jd_k = fdcomm.DLprecoders{jj,k};
    U_jd_k = fdcomm.DL_WMMSE_RX{jj,k};
    W_jd_k = fdcomm.DL_weights{jj,k};
    napla_PiB_R_DL_sum = napla_PiB_R_DL_sum-mu_gd_k*2*H_i_j'/R_in_d_j_k*H_Bj*P_jd_k/...
        (eye(d_DL(jj))+P_jd_k'*H_Bj'/R_in_d_j_k*H_Bj*P_jd_k)*...
        P_jd_k'*H_Bj'/R_in_d_j_k*H_i_j*tilde_P_iu_k;
    alpha_j = alpha_DL(jj);
    inverse_DL = inverse_DL+alpha_j*H_i_j'*U_jd_k'*W_jd_k*U_jd_k*H_i_j;
end
N_i = fdcomm.UL_UE_Ant(ii);
inverse = inverse_UL+inverse_DL + lambda_iu_k*eye(N_i);
%% Derivative of Xi_r
Nr = radar.RX;
eta_i = fdcomm.ULchannelgains(ii);
temp1 = 0;
n = radar.CUT_Idx;
d_iB = fdcomm.ULsymbols{ii};
d_iB_k_n = d_iB(:,n+1,k);
K = radar.codelength;
for nr = 1:Nr
    U_r_nr = radar.WMMSE_RX{nr,1};
    u_r_nr_k = U_r_nr(:,k);
    W_r_nr = radar.WMMSE_weights{nr,1};
    for m = 1:K
        P_iB_m = fdcomm.ULprecoders{ii,m};
        d_iB_m_n = d_iB(:,n+1,m);
        u_r_nr_m = U_r_nr(:,m);
        temp1 = temp1 + P_iB_m*d_iB_m_n*u_r_nr_m'*W_r_nr*u_r_nr_k*d_iB_k_n';
    end
end
napla_P_iB_k_Xir = 2*eta_i^2*temp1;
%% Calculate PiB
P_iu_k = inverse\(alpha_iu*H_iB'*U_iu_k'*W_iu_k+...
    (napla_PiB_R_DL_sum+napla_PiB_R_UL_sum-napla_P_iB_k_Xir)/2);
fdcomm.ULprecoders{ii,k} = P_iu_k;
end
