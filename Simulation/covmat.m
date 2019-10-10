function [cov] = covmat(fdcomm,radar,radar_comm)
%
%
%% radar only covariances
Mr                          = radar.TX; % Number of radar TX antennas
Nr                          = radar.RX; % Number of radar RX antennas
K                           = radar.codelength; % The length of the radar code; or the number of PRIs
N                           = radar.PRI_num;% Number of range cells
n                           = radar.CUT_Idx; % CUT index
%H_rtr                       = radar.channel;
Qr                          = radar.doppler; %Mr*K*Nr radar temporal steering vector
R_C                         = radar.cluttercov;
eta_t                       = radar.channelgain;
A = radar.codematrix;
sigma_0 = radar.noisepower;

%% FD comm related covariance matrices
Mc                          = fdcomm.BSTX;% Number of BS TX antennas
Nc                          = fdcomm.BSRX;% Number of BS RX antennas
I                           = fdcomm.UL_num; % Number of UL UEs
J                           = fdcomm.DL_num; % Number of DL UEs
N_UL                        = fdcomm.UL_UE_Ant; % Number of the UL UE antennas
N_DL                        = fdcomm.DL_UE_Ant; % Number of the DL UE antennas
%Htr                         = zeros(Mr,Nr);
d_DL                        = fdcomm.DLstream_num;% number of stream for each DL UE d_DL(i) <= N_DL(i)   
d_UL                        = fdcomm.ULstream_num;
P_BJ                        = fdcomm.DLprecoders; % cell(J,K*N)
P_IB                        = fdcomm.ULprecoders; % cell(I,k*N)
H_UL_DL                     = fdcomm.ULDLchannels;
D_DL                        = fdcomm.DLsymbols;
D_UL                        = fdcomm.ULsymbols;
H_DL                        = fdcomm.DLchannels;
H_UL                        = fdcomm.ULchannels;
H_BB                        = fdcomm.BBchannel;
eta_BB                      = fdcomm.BBchannelgain;
R_IB = cell(I,K);   % All the UL covariance matrices
R_BJ = cell(J,K);   % All the DL covarianc matrices
R_MUI_DL = cell(J,K*N);
R_MUI_UL = cell(I,K*N);
R_sum_DL = cell(K*N,1);
R_sum_UL = cell(K*N,1);
R_BB = cell(K*N,1);
for k = 1:K
    R_sum = zeros(Nc,Nc);
    for ii = 1 : I
        H_iB = H_UL{ii}; %load the UL channel matrix
        PiB_k = P_IB{ii,k};
        R_IB{ii,k} = H_iB*(PiB_k*PiB_k')*H_iB';
        R_sum = H_iB*(PiB_k*PiB_k')*H_iB'+R_sum;
    end
    R_sum_UL{k} = R_sum;
end
for k = 1:K
    R_sum_k = R_sum_UL{k};
    for ii = 1:I
        R_MUI_UL{ii,k} = R_sum_k - R_IB{ii,k};
    end
end
for k = 1 : K
    for ii = 1:J
        HBj = H_DL{ii}; %load the UL channel matrix
        PBj_k = P_BJ{ii,k};
        R_BJ{ii,k} = HBj*(PBj_k*PBj_k')*HBj'; 
    end
end
for k = 1 : K
    for jj = 1:J
        HBj = H_DL{jj}; %load the UL channel matrix
        jj_prime = [1:jj-1 jj+1:J];
        Nj = N_DL(jj);
        R_j_MUI = zeros(Nj,Nj);
        for jjj = 1:length(jj_prime)
            j_mui = jj_prime(jjj);
            Pj_mui = P_BJ{j_mui,k};
            R_j_MUI = HBj*(Pj_mui*Pj_mui')*HBj'+R_j_MUI;
        end
        R_MUI_DL{jj,k} = R_j_MUI;
    end
end
cov.MUI_DL = R_MUI_DL;
cov.MUI_UL = R_MUI_UL;
for k = 1 : K
    for jj = 1:J
        R_bb = zeros(Nc,Nc);
        PBj_k = P_BJ{jj,k};
        R_bb = PBj_k*PBj_k'+R_bb;
    end
    R_BB{k,1} = H_BB*R_bb*H_BB';
end
cov.B2B = R_BB;
% UL to DL
R_ULDL = cell(J,I);
for jj = 1:J
    for kk = 1:K
        R_ULDL_temp = 0; 
        for ii = 1:I
            Hij = H_UL_DL{ii,jj}; %load the UL channel matrix
            PiB_k = P_IB{ii,k};
            R_ULDL_temp = Hij*(PiB_k*PiB_k')*Hij'+R_ULDL_temp;
        end
        R_ULDL{jj,kk} = R_ULDL_temp;
    end
end
cov.UL2DL = R_ULDL;

%% Coexistence related covariance matrices
% BS-multipath/target-radar
Q_Btr = radar_comm.Btrdoppler;
%H_Bmr = radar_comm.Bmrchannels;
eta_B = radar_comm.Bmrchannelgains;
Q_Bmr = radar_comm.Bmrdoppler;
R_Bmr = zeros(K,K,Nr);
V_Btr = zeros(K,Mc); % DL signal vectors target
V_Bmr = zeros(K,Mc); % DL signal vectors multipath 
for k = 1:K
    v_Bt_k = 0;
    v_Bm_k = 0;
    for jj = 1:J
        P_Bj_k    = P_BJ{jj,k};
        v_Bt_k  = P_Bj_k*D_DL{jj}(:,1,k) + v_Bt_k;
        v_Bm_k  = P_Bj_k*D_DL{jj}(:,n+1,k) + v_Bm_k;
    end
    V_Bmr(k,:) = v_Bm_k.';
    V_Btr(k,:) = v_Bt_k.';
end
f_Bm_Nr = radar_comm.f_Bm_Nr;
for nr = 1:Nr
    f_Bmnr = f_Bm_Nr(nr);
    for m = 1:K
        for l = 1:K
            v_Bm_m = V_Bmr(m,:);
            v_Bm_l = V_Bmr(l,:);
            R_Bmr(m,l,nr) = eta_B^2*exp(1i*2*pi*(m-l)*f_Bmnr)*conj(v_Bm_l)*v_Bm_m.';
        end
    end
end
% covariance matrix R_t_r 
R_tr = zeros(K,K,Nr); %Covariance for the radar
S_rtr = zeros(K,Mr,Nr); % signal matrix of the radar 
S_Btr = zeros(K,Mc,Nr);
S_tr = zeros(K,Mr,Nr);
for nr = 1:Nr
    f_Bt_nr = radar_comm.f_Bt_Nr(nr);
    for m = 1:K
        q_rnr_m = Qr(:,m,nr);
        Q_rnr_m = diag(q_rnr_m);
        v_Bt_m  = V_Btr(m,:);
        a_m = A(m,:);
        S_rtr(m,:,nr) = Q_rnr_m * a_m.';
        q_Btnr_m= exp(1i*2*pi*(m-1)*f_Bt_nr);
        S_Btr(m,:,nr) = q_Btnr_m*v_Bt_m;
        for l = 1:K
            q_Btnr_l= exp(1i*2*pi*(l-1)*f_Bt_nr);
            q_rnr_l = Qr(:,l,nr);
            Q_rnr_l = diag(q_rnr_l);
            v_Bt_l  = V_Btr(l,:);
            a_l = A(l,:);
            R_tr(m,l,nr) = eta_t(nr)^2*(conj(a_l)*Q_rnr_l'+conj(q_Btnr_l)*conj(v_Bt_l))*(Q_rnr_m*a_m.'+q_Btnr_m*v_Bt_m.');
        end
    end
    S_tr(:,:,nr) = S_rtr(:,:,nr) + S_Btr(:,:,nr);
end
cov.S_tr = S_tr;
cov.S_rtr = S_rtr;
cov.S_Btr = S_Btr;
cov.target2radar = R_tr;
% UL UEs to radar RX
Q_Ir = radar_comm.ULrdoppler;
%H_UL_r = radar_comm.UL2radarchannels;
eta_U = radar_comm.UL2radarchannelgains;
V_I = cell(I,1);
for ii = 1:I
    N_i = N_UL(ii);
    v_i = zeros(K,N_i);
    for k = 1:K
        P_iB_k = P_IB{ii,k};
        d_iB_k = D_UL{ii}(:,1,k);
        v_i(k,:) = d_iB_k.'*P_iB_k.';
    end
    V_I{ii} = v_i;
end
R_Ur = zeros(K,K,Nr);
for nr = 1:Nr
    for m = 1:K
        for l = 1:K
            R_U_nr_temp = 0;
            for ii = 1:I
                v_i = V_I{ii};
                v_i_m = v_i(m,:);
                v_i_l = v_i(l,:);
                f_i_nr = log(Q_Ir(2,nr,ii))/(2*pi*1i);
                R_U_nr_temp = eta_U(ii,nr)^2*exp(1i*2*pi*(m-l)*f_i_nr)*conj(v_i_l)*v_i_m.'...
                   +R_U_nr_temp;
            end
            R_Ur(m,l,nr) = R_U_nr_temp;
        end
    end
end
cov.UL2radar = R_Ur;
% radar to BS
R_rB = zeros(Mc,Mc,K);
H_r_BS = radar_comm.radar2BSchannels;
for k = 1:K 
    ak = A(k,:).';
    R_rB(:,:,k) = (H_r_BS*ak)*(ak'*H_r_BS');
end
cov.radar2BS = R_rB;
% radar to DL UEs
R_rJ    = cell(J,1);
H_r_DL  = radar_comm.radar2DLchannnels;
for jj = 1:J
    N_j = N_DL(jj);
    R_rj = zeros(N_j,N_j);
    H_rj = H_r_DL{jj};
    for k = 1:K
        ak = A(k,:).';
        R_rj(:,:,k) = H_rj*ak*(ak'*H_rj');
    end
    R_rJ{jj,1} = R_rj;
end
cov.radar2DL = R_rJ;
% radar interference matrices
R_Zr = zeros(K,K,Nr);
R_in_r = zeros(K,K,Nr);
sigma_n = radar.noisepower;
for nr = 1:Nr
    R_Zr(:,:,nr) = sigma_n*eye(K,K);
    R_in_r(:,:,nr) = R_Zr(:,:,nr)+R_Ur(:,:,nr)+R_Bmr(:,:,nr)+R_C(:,:,nr);
end
cov.noise = R_Zr;
cov.inr = R_in_r;
% UL interference matrices
R_in_UL = cell(I,K);
R_total_UL = cell(I,K);
R_ZB = sigma_0*eye(Nc,Nc);
for ii = 1:I
    for k = 1:K
        R_MUI_i_k = R_MUI_UL{ii,k};
        R_BB_k = R_BB{k,1};
        R_rB_k = R_rB(:,:,k);
        R_iB_k = R_IB{ii,k};
        R_in_UL{ii,k} = R_MUI_i_k+ R_BB_k + R_rB_k + R_ZB;
        R_total_UL{ii,k} = R_in_UL{ii,k} + R_iB_k;
    end
end
cov.in_UL = R_in_UL;
cov.total_UL = R_total_UL;
% DL interference matrices
R_in_DL = cell(J,K);
R_total_DL = cell(J,K);
for jj = 1:J
    R_rj = R_rJ{jj,1};
    N_j = N_DL(jj);
    for k = 1:K
        R_MUI_j_k = R_MUI_DL{jj,k};
        R_UL_j_k  = R_ULDL{jj,k};
        R_rj_k = R_rj(:,:,k);
        R_in_DL{jj,k} = R_MUI_j_k + R_UL_j_k + R_rj_k + sigma_n*eye(N_j,N_j);
        R_dj_k = R_BJ{jj,k};
        R_total_DL{jj,k} = R_dj_k + R_MUI_j_k + R_UL_j_k + R_rj_k + sigma_n*eye(N_j,N_j);
    end
end
cov.in_DL = R_in_DL;
cov.total_DL = R_total_DL;

end

