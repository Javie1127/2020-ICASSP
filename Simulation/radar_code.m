function [radar] = radar_code(fdcomm, radar, radar_comm, cov,k)
% radar code matrix update
%
I = fdcomm.UL_num;
J = fdcomm.DL_num;
Mr = radar.TX;
K = radar.codelength;
Nr = radar.RX;
n = radar.CUT_Idx;
%eta_B = radar_comm.Bmrchannelgains;
eta_T = radar.channelgain;
H_rB = radar_comm.radar2BSchannels;
%% matrix D
lambda_r_k = radar.lambda(k);
D = 0;
for ii = 1 : I
    U_iu_k = fdcomm.UL_WMMSE_RX{ii,k};
    W_iu_k = fdcomm.UL_weights{ii,k};
    D = D +  H_rB'*U_iu_k'*W_iu_k'*U_iu_k*H_rB;
end
for jj = 1:J
    U_jd_k = fdcomm.DL_WMMSE_RX{jj,k};
    W_jd_k = fdcomm.DL_weights{jj,k};
    H_rj = radar_comm.radar2DLchannnels{jj};
    D = D + H_rj'*U_jd_k'*W_jd_k*U_jd_k*H_rj;
end
D = D + lambda_r_k*eye(Mr);
%% Matrix E 
E = 0;
for nr = 1:Nr
    Urnr = radar.WMMSE_RX{nr,1};
    urnr_k = Urnr(:,k);
    Wrnr = radar.WMMSE_weights{nr};
    E = urnr_k'*Wrnr*urnr_k + E;
end
%% Matrix F 
d_UL = fdcomm.ULstream_num; % number of data streams of the UL UEs
d_DL = fdcomm.DLstream_num; % number of data streams of the DL UEs
napla_ak_R_UL_sum = 0;
tilde_ar_k = radar.codematrix(k,:).';
for ii = 1 : I
    mu_iu_k = fdcomm.mu_UL(ii,k);
    R_in_iu_k = cov.in_UL{ii,k};
    H_iB = fdcomm.ULchannels{ii,1};
    P_iu_k = fdcomm.ULprecoders{ii,k};
    napla_ak_R_UL_sum = napla_ak_R_UL_sum - mu_iu_k*2*H_rB'/R_in_iu_k*H_iB*P_iu_k/...
        (eye(d_UL(ii))+P_iu_k'*H_iB'/R_in_iu_k*H_iB*P_iu_k)*...
        P_iu_k'*H_iB'/R_in_iu_k*H_rB*tilde_ar_k;
end
napla_ak_R_DL_sum = 0;
for jj = 1:J
    mu_jd_k = fdcomm.mu_DL(jj,k);
    H_rj = radar_comm.radar2DLchannnels{jj};
    R_in_jd_k = cov.in_DL{jj,k};
    H_Bj = fdcomm.DLchannels{jj,1};
    P_jd_k = fdcomm.DLprecoders{jj,k};
    napla_ak_R_DL_sum = napla_ak_R_DL_sum-mu_jd_k*2*H_rj'/R_in_jd_k*H_Bj*P_jd_k/...
        (eye(d_DL(jj))+P_jd_k'*H_Bj'/R_in_jd_k*H_Bj*P_jd_k)*...
        P_jd_k'*H_Bj'/R_in_jd_k*H_rj*tilde_ar_k;
end
Qr = radar.doppler;
term = 0;
f_Bt_Nr = radar_comm.f_Bt_Nr;
for nr = 1:Nr
    eta_t = eta_T(nr);
    Urnr = radar.WMMSE_RX{nr,1};
    urnr_k = Urnr(:,k);
    Wrnr = radar.WMMSE_weights{nr,1};
    q_Bt_nr = radar_comm.Btrdoppler(:,nr);
    f_Bt_nr = f_Bt_Nr(nr);
    q_r_nr_k = Qr(:,k,nr);
    Q_rnr_k = diag(q_r_nr_k);
    temp_1 = 0;
    for g = 1:J
        d_gd = fdcomm.DLsymbols{g,1};
        d_gd_k_1 = d_gd(:,1,k);
        P_gd_k = fdcomm.DLprecoders{g,k};
        temp_1 = P_gd_k*d_gd_k_1+temp_1;
    end
    term_1 = Q_rnr_k'-cos(2*pi*(k-1)*f_Bt_nr)*Q_rnr_k'*temp_1*urnr_k';
    term_1 = eta_t^2*term_1;
    term_2 = 0;
    for m = 1:K
        if m ~= k
            urnr_m = Urnr(:,m);
            q_r_nr_m = Qr(:,m,nr);
            Q_r_nr_m = diag(q_r_nr_m);
            ar_m = radar.codematrix(m,:).';
            temp_2 = 0;
            for g = 1:J
                d_gd = fdcomm.DLsymbols{g,1};
                d_gd_m_n_1 = d_gd(:,n+1,m);
                P_gd_m = fdcomm.DLprecoders{g,m};
                temp_2 = P_gd_m*d_gd_m_n_1 + temp_2;
            end
            temp_2 = q_Bt_nr(m)*temp_2+ Q_r_nr_m*ar_m;
            term_2 = eta_t^2*Q_rnr_k'*temp_2*urnr_m'+term_2;
        term = (term_1 + term_2)*Wrnr*urnr_k + term;
        end
    end
F = napla_ak_R_DL_sum+napla_ak_R_UL_sum+term;
a_k = sylvester(D,E,F);
radar.codematrix(k,:) = a_k.'; 
end
end

