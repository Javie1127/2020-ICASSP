function [fdcomm,radar,radar_comm] = parameters_icassp
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Array Parameters
Mr                          = 8; % Number of radar TX antennas
Mc                          = Mr;% Number of BS TX antennas
Nr                          = 6; % Number of radar RX antennas
Nc                          = Mc;% Number of BS RX antennas
I                           = 5; % Number of UL UEs
J                           = 4; % Number of DL UEs
K                           = 7; % The length of the radar code; or the number of PRIs
N_UL                        = 4*ones(I,1); % Number of the UL UE antennas
N_DL                        = 3*ones(J,1); % Number of the DL UE antennas
%Htr                         = zeros(Mr,Nr);
d_DL                        = 3*ones(J,1);% number of stream for each DL UE d_DL(i) <= N_DL(i) 
d_UL                        = 4*ones(I,1); 
N                           = 10;% Number of range cells
n                           = 4; % CUT index
%% Power 
sigma_0                     = 0.1; % noise power
SNR_rtr                     = 0.1*ones(Mr,Nr); % db radar TX signal to noise ratio at the radar RX
SNR_BS_DL                   = 1*ones(J,1);
SNR_UL_BS                   = 2*ones(I,1);
P_r                         = 2*ones(K,1);
P_DL                        = 0;
P_UL                        = zeros(I,1);
%% radar-only paramenters
H_rtr                       = zeros(Mr,Nr);
eta_t                       = ones(Nr,1);
%eta_t                       = 1;
for mr = 1 : Mr
    P_r(mr) = sum(db2pow(SNR_rtr(mr,:))*sigma_0);
    for nr = 1:Nr
        H_rtr (mr,nr) = sqrt(eta_t(nr)/2)*(randn(1)+1i*randn(1));
    end
end
Qr = zeros(Mr,K,Nr);
for nr = 1:Nr
    for k = 1 : K
        % the model of f_tnr is the same as ...
        % "An Information Theoretic Approach to Robust Constrained Code Design for MIMO Radars"
        f_mr_nr              = -1*rand(Mr,1)+0.5;  %[-0.5,0.5]
        f_mr_t_nr           = 0.15*f_mr_nr+0.25; % Normalized Doppler frequency
        Qr(:,k,nr)          = exp(1j*2*pi*(k-1).*f_mr_t_nr); % Doppler Domain Steering vector 
    end
end
% Clutter 
CNR                             = 0*ones(Nr,1);
R_C                             = zeros(K,K,nr); 
rho = 0.5;
% Method 1 C is a gaussian rv (0,sigma_c);
for nr = 1:Nr
    sigma_c_nr = db2pow(CNR(nr))*sigma_0;
    for ii = 1 : K 
        for jj = 1:K
            R_C(ii,jj,nr) = sigma_c_nr^2*rho^(ii-jj);
        end
    end
end
radar.channel               = H_rtr;
radar.doppler                = Qr; %radar temporal steering vector
radar.cluttercov            = R_C;
radar.channelgain           = eta_t;
radar.power                 = P_r;
radar.TX = Mr;
radar.RX = Nr;
radar.codelength = K;
radar.PRI_num = N;
radar.CUT_Idx = n;
radar.noisepower = sigma_0;
%%% WMMSE receiver and weights 
radar.WMMSE_RX = cell(Nr,1);
radar.WMMSE_weights = cell(Nr,1);
radar.MMSE = cell(Nr,1);
radar.MMSE_nop = cell(Nr,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);
radar.Xi_r = zeros(Nr,1);
radar.lambda = ones(K,1);
%-----------------------------------------------------------------%
%% FD Comm
% UL channels
eta_UL                      = ones(I,1);  % Channel gains of the UL
H_UL                        = cell(I,1);  %
for ii = 1:I
    P_UL(ii)                = db2pow(SNR_UL_BS(ii))*sigma_0; 
    H_UL{ii,1}                =...
        sqrt(eta_UL(ii)/2)*(randn(Nc,N_UL(ii))+1i*randn(Nc,N_UL(ii))); 
end
% DL Channels
H_DL                        = cell(J,1); % J DL Channels
eta_DL                      = ones(J,1); %
for jj = 1:J
    sigma_BS_j              = db2pow(SNR_BS_DL(jj))*sigma_0;
        P_DL                = P_DL + sigma_BS_j;
    H_DL{jj,1}                = sqrt(eta_DL(jj)/2/N_DL(jj))*(randn(N_DL(jj),Mc)+1i*randn(N_DL(jj),Mc)); 
end
% communication qpsk symbols
mod =   comm.QPSKModulator;
%refC =  constellation(mod);
D_DL = cell(J,1);
D_UL = cell(I,1);
for jj = 1:J
    d_j = d_DL(jj);
    frames = zeros(d_j, N, K);
    for kk = 1:K
        for nn = 1:N
            v = randn(d_j,2);
            v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));
            frames(:,nn,kk) = complex(v(:,1),v(:,2));
        end
    end
    D_DL{jj,1} = frames;
end
for ii = 1:I
    d_i = d_UL(ii);
    frames = zeros(d_i, N, K);
    for kk = 1:K
        for nn = 1:N
            v = randn(d_i, 2);
            v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));
            frames(:,nn,kk) = complex(v(:,1),v(:,2));
        end
    end
    D_UL{ii,1} = frames;
end
% UL UE - DL UE channel
eta_UL_DL                   = ones(I,J);  
H_UL_DL                     = cell(I,J);
for ii = 1:I
    for jj = 1:J
        H_UL_DL{ii,jj}      = sqrt(eta_UL_DL(ii,jj)/2)*(randn(N_DL(jj),N_UL(ii))+1i*randn(N_DL(jj),N_UL(ii)));
    end
end
% BS-to-BS channel
eta_BB                      = 1;  % self-interference channel gain 
H_BB                        = sqrt(eta_BB/2)*(randn(Nc,Mc)+1i*randn(Nc,Mc));
fdcomm.DLsymbols        = D_DL;
fdcomm.ULsymbols        = D_UL;
fdcomm.DLchannels       = H_DL;
fdcomm.ULchannels       = H_UL;
fdcomm.DLchannelgains   = eta_DL;
fdcomm.ULchannelgains   = eta_UL;
fdcomm.DLpower          = P_DL;
fdcomm.ULpower          = P_UL;
fdcomm.ULDLchannels     = H_UL_DL;
fdcomm.ULDLchannelgains = eta_UL_DL;
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.UL_UE_Ant = N_UL;
fdcomm.DL_UE_Ant = N_DL;
fdcomm.DLstream_num = d_DL;
fdcomm.ULstream_num = d_UL;
fdcomm.BBchannel = H_BB;
fdcomm.BBchannelgain = eta_BB;
fdcomm.UL_WMMSE_RX = cell(I,K);
fdcomm.DL_WMMSE_RX = cell(J,K);
fdcomm.UL_weights = cell(I,K);
fdcomm.DL_weights = cell(J,K);
fdcomm.UL_MMSE = cell(I,K);
fdcomm.DL_MMSE = cell(J,K);
fdcomm.UL_MMSE_nop = cell(I,K);
fdcomm.DL_MMSE_nop = cell(J,K);
fdcomm.alpha_UL = 1/(I+J+Nr)*ones(I,1);
fdcomm.alpha_DL = 1/(I+J+Nr)*ones(J,1);


%% Coexistence
% BS- target- radar RX
%SNR_Btr = 10*ones(Nr,1); % db BS signal to noise ratio at the nrth radar RX
%H_Btr   = zeros(Mc,Nr);
%eta_Btr = ones(Nr,1);
% eta_Bt   = 1; 
% for nr = 1 : Nr
%     H_Btr(:,nr) = sqrt(eta_Bt/2).*(randn(Mc,1)+1i*randn(Mc,1));
% end
% Q_Btr = zeros(K,Nr);
% for nr = 1:Nr
%     f_B_t_nr            = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
%     Q_Btr(:,nr)            = exp(1j*2*pi*f_B_t_nr.*kk'); % Doppler Domain Steering vector 
% end
% BS - multipath - radar RX
kk = 0:K-1;
eta_B = 1;
f_Bm_Nr = zeros(Nr,1);
Q_Btr = zeros(K,Nr);
Q_Bmr = zeros(K,Nr);
Q_Ir  = zeros(K,Nr,I);
H_Bmr = zeros(Mc,Nr);
f_Bt_Nr = zeros(Nr,1);
for nr = 1:Nr
    f_Bt_nr                  = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
    f_Bm_nr                  = 0.15*(-1*rand(1)+0.5)+0.25;
    Q_Btr(:,nr)              = exp(1j*2*pi*f_Bt_nr.*kk'); % Doppler Domain Steering vector 
    Q_Bmr(:,nr)              = exp(1j*2*pi*f_Bm_nr.*kk');
    H_Bmr(:,nr)              = sqrt(eta_B/2).*(randn(Mc,1)+1i*randn(Mc,1));
    f_Bt_Nr(nr)              = f_Bt_nr;
    f_Bm_Nr(nr)              = f_Bm_nr;
end
for ii = 1:I
    for nr = 1:Nr
        f_i_nr = 0.15*(-1*rand(1)+0.5)+0.25;
        Q_Ir(:,nr,ii) = exp(1j*2*pi*f_i_nr.*kk');
    end
end
% UL UEs - radar RX
INR_UL_r = ones(I,Nr);
H_UL_r   = cell(I,Nr);
for ii = 1:I
    for nr = 1:Nr
        sigma_i_nr = db2pow(INR_UL_r(ii,nr));
        H_UL_r{ii,nr} = sqrt(sigma_i_nr/2)*(randn(N_UL(ii),1)+1i*randn(N_UL(ii),1));
    end
end
% radar - DL UE
eta_r_DL = ones(Mr,J);
%INR_r_DL = ones(Mr,J);
H_r_DL = cell(J,1);
mu_DL = 0.2*ones(Mr,J);
kappa_DL = 0.4;
for jj = 1:J
    N_j = N_DL(jj);
    H_r_j = zeros(N_j,Mr);
    for mr = 1 : Mr
        %sigma_mr_j = db2pow(INR_r_DL(mr,jj));  
        H_r_j(:,mr) = sqrt(eta_r_DL(mr,jj)^2/(kappa_DL+1)/2).*(randn(N_j,1)+1i*randn(N_j,1))...
            +sqrt(kappa_DL/(kappa_DL+1))*mu_DL(mr,jj).*ones(N_j,1);
    end
    H_r_DL{jj} = H_r_j;
end
% radar - BS Rx
%eta_r_BS = ones(Mr,1);
%INR_r_BS = ones(Mr,1);
H_r_BS = zeros(Nc,Mr);
mu_Mr_BS = 0.5*ones(Mr,1);
kappa_BS = 0.5;
eta_Mr = 1*ones(Mr,1);
for mr = 1 : Mr
    %sigma_mr_BS = db2pow(INR_r_DL(mr));
    H_r_BS(:,mr) = sqrt(eta_Mr(mr)/(kappa_BS+1)/2).*(randn(Nc,1)+1i*randn(Nc,1))+sqrt(kappa_BS/(kappa_BS+1))*mu_Mr_BS(mr)*ones(Nc,1);
end
% UL UEs - Radar Rxs
eta_U = ones(I,Nr);
H_UL_r = cell(I,Nr);
for ii = 1:I
    for nr = 1:Nr
        H_UL_r{ii,nr} = sqrt(eta_U(ii,nr)/2).*(randn(N_UL(ii),1)+1i*randn(N_UL(ii),1));
    end
end
radar_comm.f_Bt_Nr = f_Bt_Nr;
radar_comm.f_Bm_Nr = f_Bm_Nr;
radar_comm.Btrdoppler = Q_Btr;
radar_comm.Bmrchannels = H_Bmr;
radar_comm.Bmrchannelgains = eta_B;
radar_comm.Bmrdoppler = Q_Bmr;
radar_comm.ULrdoppler = Q_Ir;
radar_comm.UL2radarchannels = H_UL_r;
radar_comm.UL2radarchannelgains = eta_U;
radar_comm.radar2BSchannels = H_r_BS;
radar_comm.radar2BSchannelgains = eta_Mr;
radar_comm.radar2BSchannelmeans = mu_Mr_BS;
radar_comm.radar2BSchannelkfactor = kappa_BS;
radar_comm.radar2DLchannnels = H_r_DL;
radar_comm.radar2DLchannnelgains = eta_r_DL;
%% Alternating optimization
radar_comm.iter1_max = 10;
radar_comm.iter2_max = 10;
end
