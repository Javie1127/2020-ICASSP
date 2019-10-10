function [fdcomm] = DL_precoders(k, jj, fdcomm, radar_comm, radar,cov)
%function [P_iB_k] = UL_precoders(ii, k,...
%   H_UL, H_i_MUI, W_DL_k, W_UL_k, U_DL_k, U_UL_k,Ur,...
%   lambda_UL, mu_DL, mu_UL, napla_PiB_R_DL_k,...,
%    napla_PiB_R_UL_k, napla_PiBk_Xi_r)
% H_DL: downlink channel matrices
% H_BB: Self-interference channel matrix
% W_DL_k: DL UE Weight Matices
% W_UL_k: UL UE Weight Matrices
% Wr : radar Weight matrices; 
% A: K*Mr
% F_rtr Mr*Nr contains the doppler frequencies of each TX-RX pair
% lambda_UL: Lagrange multiplier vector for the UL power constraints
% mu_DL: K*J Lagrange multiplier vector for the DL rate constraints
% mu_UL: K*I Lagrange multiplier vector for the UL rate constraints
% napla_PBj_R_DL_k: Mc*dj*J
% napla_PBj_R_UL_k: Mc*dj*I
I = fdcomm.UL_num;
J = fdcomm.DL_num;
Nr = radar.RX;
K = radar.codelength;
d_DL = fdcomm.DLstream_num; % number of data streams of the DL UEs
d_UL = fdcomm.ULstream_num; % number of data streams of the UL UE
%% Matrix A 
A_DL = 0;
for g = 1:J
    H_Bg = fdcomm.DLchannels{g,1};
    U_Bg_k = fdcomm.DL_WMMSE_RX{g,k};
    W_Bg_k = fdcomm.DL_weights{g,k};
    A_DL = H_Bg'*U_Bg_k'*W_Bg_k*U_Bg_k*H_Bg+A_DL;
end
A_UL =0;
H_BB = fdcomm.BBchannel;
Mc = fdcomm.BSTX;
for ii = 1:I
    U_iu_k = fdcomm.UL_WMMSE_RX{ii,k};
    W_iu_k = fdcomm.UL_weights{ii,k};
    A_UL = A_UL + H_BB'*U_iu_k'*W_iu_k*U_iu_k*H_BB;
end
lambda_k_d = fdcomm.lambda_DL(k);
A = A_DL + A_UL +lambda_k_d*eye(Mc);
%% Matrix B
d_jd = fdcomm.DLsymbols{jj,1};
n = radar.CUT_Idx;
d_jd_k_n = d_jd(:,n+1,k);
d_jd_k_1 = d_jd(:,1,k);
eta_B = radar_comm.Bmrchannelgains;
eta_t = radar.channelgain;
d_j = fdcomm.DLstream_num(jj);
B = zeros(d_j,d_j);
for nr = 1: Nr
    Urnr = radar.WMMSE_RX{nr,1};
    urnr_k = Urnr(:,k);
    Wrnr = radar.WMMSE_weights{nr};
    B_temp = eta_B^2*d_jd_k_n* urnr_k'*Wrnr*urnr_k*d_jd_k_n'+...
       eta_t(nr)^2*d_jd_k_1*urnr_k'*Wrnr*urnr_k*d_jd_k_1';
   B = B_temp+B;
end
%% Matrix C
% derivative of downlink rate
napla_P_jd_k_R_DL_k = 0;
tilde_P_jd_k = fdcomm.DLprecoders{jj,k};
H_Bj = fdcomm.DLchannels{jj,1};
for g = 1:J
    R_in_gd_k = cov.in_DL{g,k};
    mu_gd_k = fdcomm.mu_DL(g,k);
    if g == jj
       napla_P_jd_k_R_DL_k = napla_P_jd_k_R_DL_k + mu_gd_k*2*H_Bj'/R_in_gd_k*H_Bj*tilde_P_jd_k/...
           (eye(d_DL(g))+tilde_P_jd_k'*H_Bj'/R_in_gd_k*H_Bj*tilde_P_jd_k);
    else
        H_Bg = fdcomm.DLchannels{g,1};
        P_gd_k = fdcomm.DLprecoders{g,k};
        napla_P_jd_k_R_DL_k = napla_P_jd_k_R_DL_k - mu_gd_k*2*H_Bg'/R_in_gd_k*H_Bg*tilde_P_jd_k/...
           (eye(d_DL(g))+P_gd_k'*H_Bg'/R_in_gd_k*H_Bg*P_gd_k)*...
           P_gd_k'*H_Bg'/R_in_gd_k*H_Bg*tilde_P_jd_k;
    end
end
% derivative sum w.r.t. the UL rates
napla_P_jd_k_R_UL_k = 0;
for q = 1:I
    mu_qu_k = fdcomm.mu_DL(g,k);
    H_iB = fdcomm.ULchannels{q,1};
    R_in_iu_k = cov.in_UL{q,k};
    P_iu_k = fdcomm.ULprecoders{q,k};
    napla_P_jd_k_R_UL_k = napla_P_jd_k_R_UL_k - mu_qu_k*2*H_BB'/R_in_iu_k*H_iB*P_iu_k/...
        (eye(d_UL(q))+P_iu_k'*H_iB'/R_in_iu_k*H_iB*P_iu_k)*...
        P_iu_k'*H_iB'/R_in_iu_k*H_BB*tilde_P_jd_k;
end
U_jd_k = fdcomm.DL_WMMSE_RX{jj,k};
W_jd_k = fdcomm.DL_weights{jj,k};
CC = H_Bj'*U_jd_k'*W_jd_k+napla_P_jd_k_R_UL_k+napla_P_jd_k_R_DL_k;
f_Bt_Nr = radar_comm.f_Bt_Nr;
Qr = radar.doppler;
ar_k = radar.codematrix(k,:).';
term = 0;
for nr = 1 : Nr
    Urnr = radar.WMMSE_RX{nr,1};
    urnr_k = Urnr(:,k);
    Wrnr = radar.WMMSE_weights{nr,1};
    q_Bt_nr = radar_comm.Btrdoppler(:,nr);
    f_Bt_nr = f_Bt_Nr(nr);
    q_r_nr_k = Qr(:,k,nr);
    Q_rnr_k = diag(q_r_nr_k);
    term_1 = eta_t(nr)^2*(q_Bt_nr(k)-cos(2*pi*(k-1)*f_Bt_nr)*Q_rnr_k*ar_k*urnr_k')*Wrnr*urnr_k*d_jd_k_1';
    term_2 = 0;
    term_3 = 0;
    for g = 1:J
        if g ~= jj
            P_gd_k = fdcomm.DLprecoders{g,k};
            d_gd = fdcomm.DLsymbols{g,1};
            d_gd_k_n_1 = d_gd(:,n+1,k);
            term_2 = P_gd_k*d_gd_k_n_1*urnr_k'*Wrnr*urnr_k*d_gd_k_n_1'+term_2;
            d_gd_k_1 = d_gd(:,1,k);
            term_3 = eta_t(nr)^2*P_gd_k*d_gd_k_1*urnr_k'*Wrnr*urnr_k*d_gd_k_1' + term_3;
        end
    end
    term_2 = eta_B^2*term_2;
    term_4 = 0;
    term_5 = 0;
    for m = 1:K
        if m ~= k
            urnr_m = Urnr(:,m);
            q_r_nr_m = Qr(:,m,nr);
            Q_r_nr_m = diag(q_r_nr_m);
            ar_m = radar.codematrix(m,:).';
            temp_1 = 0;
            temp_2 = 0;
            for g = 1:J
                d_gd = fdcomm.DLsymbols{g,1};
                d_gd_m_1 = d_gd(:,1,m);
                d_gd_m_n_1 = d_gd(:,n+1,m);
                P_gd_m = fdcomm.DLprecoders{g,m};
                temp_1 = P_gd_m*d_gd_m_1+temp_1;
                temp_2 = P_gd_m*d_gd_m_n_1 + temp_2;
            end
            term_4 = (Q_r_nr_m*ar_m+ temp_1)*urnr_m'*Wrnr*urnr_k*d_jd_k_1'+term_4;
            term_5 = temp_2*urnr_m'*Wrnr*urnr_k*d_jd_k_n'+term_5;
        end
    end
    term = eta_t(nr)^2*term_1-eta_B^2*term_2-eta_t(nr)^2*term_3-eta_t(nr)^2*term_4 -eta_B^2*term_5 + term;
end
C = term + CC;
P_jd_k = sylvester(A,B,C);
fdcomm.DLprecoders{jj,k} = P_jd_k;
end
