function [fdcomm] = napla_Rate(jj, ii,k,fdcomm,radar_comm,cov)

% function [napla_PiB_R_UL_k,napla_PiB_R_DL_k,...
%     napla_PBj_R_DL_k,napla_PBj_R_UL_k,...
%  napla_a_R_UL_k,napla_a_R_DL_k]...
%     = napla_Rate(ii,jj,tilde_P_iB_k,tilde_P_Bj_k,tilde_ak,...
%     H_DL, H_UL, H_i_MUI,H_BB, H_r_DL, H_rB, ...
%     R_in_DL_k,R_in_UL_k,P_DL_k,P_UL_k)
% napla_Rate returns 
% H_DL:
% R_in_DL_k: the interference-plus-noise covariance matrices...
%            for the DL UEs during the kth frame
% R_in_UL_k: the interference-plus-noise covariance matrices...
%            for the UL UEs during the kth frame
% P_DL J*K cell array contains the current DL precoding matrices
% P_UL I*K cell array contains the current UL precoding matrices
% I = size(H_UL,3); % Number of the UL UEs
% J = size(H_DL,3); % Number of the DL UEs
%%
%--------Initialization------------------
napla_PiB_R_UL_k = cell(I,1);
napla_PiB_R_DL_k = cell(J,1);
napla_PBj_R_DL_k = cell(J,1);
napla_PBj_R_UL_k = cell(I,1);
H_iB = fdcomm.ULchannels(ii,1);
H_BB = fdcomm.BBchannel;
%% 
for g = 1:fdcomm.DL_num
    R_in_dg_k = cov.in_DL{g,k};
    H_Bg = fdcomm.ULchannels{g,k};
    P_Bg_k = fdcomm.DLprecoder{g,k};
    d_dg = size(P_Bg_k,2); % number of data streams of the gth DL UE
    H_rj = radar_comm.radar2DLchannels(g,1);
    if g == jj
        napla_PBj_R_DL_k{g} = H_Bg'/R_in_dg_k*H_Bg*tilde_P_Bj_k/...
           (eye(d_dg)+tilde_P_Bj_k'*H_Bg'/R_in_dg_k*H_Bg*tilde_P_Bj_k);
    else
        H_Bg = H_DL{g};
        napla_PBj_R_DL_k{g} = -H_Bg'/R_in_dg_k*H_Bg*P_Bg_k/...
            (eye(d_dg)+P_Bg_k'*H_Bg'/R_in_dg_k*H_Bg*P_Bg_k)*...
            P_Bg_k'*H_Bg'/R_in_dg_k*H_Bg*tilde_P_Bj_k;
    end
    Hig = H_i_MUI{g};
    napla_PiB_R_DL_k{g} = -Hig'/R_in_dg_k*HBj*P_Bg_k/...
        (eye(d_dg)+PBj'*H_Bg/R_in_dj_k*H_Bg*P_Bg_k)*...
        P_Bg_k'*H_Bg'/R_in_dg_k*Hig*tilde_P_iB_k;
        napla_a_R_DL_k = -H_rj'/R_in_dg_k*H_Bg*P_Bg_k/...
        (eye(d_dg)+PBj'*H_Bg/R_in_dj_k*H_Bg*P_Bg_k)*...
        P_Bg_k'*H_Bg'/R_in_dg_k*H_rj*tilde_ak;
end
for q = 1:I
    R_in_uq_k = R_in_UL_k{q};
    H_qB = H_UL{q};
    P_qB_k = P_UL_k{q};
    d_uq = size(P_qB_k,2); % number of data streams of the qth UL UE
    if q == ii
        napla_PiB_R_UL_k{q} = H_qB'/R_in_uq_k*H_qB*tilde_P_Bj_k/...
           (eye(d_ug)+tilde_P_Bj_k'*H_qB'/R_in_dg_k*H_Bg*tilde_P_iB_k);
    else
        napla_PiB_R_UL_k{q} = -H_qB'/R_in_uq_k*H_qB*P_qB_k/...
           (eye(d_ug)+tilde_P_Bj_k'*H_qB'/R_in_dg_k*H_Bg*P_qB_k)*...
           P_qB_k*H_qB/R_in_uq_k*H_iB*tilde_P_iB_k;
    end
    napla_PBj_R_UL_k{q} = -H_BB'/R_in_uq_k*H_qB*P_qB_k/...
        (eye(d_uq)+P_qB_k'*H_qB'/R_in_uq_k*H_qB*P_qB_k)*...
        P_qB_k'*H_qB'/R_in_uq_k*H_BB*tilde_P_Bj_k;
    napla_a_R_UL_k = -H_rB'/R_in_uq_k*H_qB*P_qB_k/...
    (eye(d_uq)+P_qB_k'*H_qB'/R_in_uq_k*H_qB*P_qB_k)*...
    P_qB_k*H_qB/R_in_uq_k*H_rB*tilde_ak;
end
end

