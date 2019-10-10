function [fdcomm] = Comm_MMSE(fdcomm, radar, cov)
%UComm_MMSE_receiver returns:
%----linear UL/DL MMSE receiver
%----Optimal Weight Matrices
K = radar.codelength;
I = fdcomm.UL_num;
J = fdcomm.DL_num;
U_UL = fdcomm.UL_WMMSE_RX;
U_DL = fdcomm.DL_WMMSE_RX;
W_UL = fdcomm.UL_weights;
W_DL = fdcomm.DL_weights;
E_UL_star = fdcomm.UL_MMSE;
E_DL_star = fdcomm.DL_MMSE;
E_UL = fdcomm.UL_MMSE_nop;
E_DL = fdcomm.DL_MMSE_nop;
for k = 1:K 
    for ii = 1 : I
        HiB = fdcomm.ULchannels{ii};
        R_UL_i_k = cov.total_UL{ii,k};
        P_iB_k = fdcomm.ULprecoders{ii,k};
        U_iB_k = P_iB_k'*HiB'/(R_UL_i_k);
        E_iB_k_star = abs(eye(fdcomm.ULstream_num(ii))-P_iB_k'*HiB'/R_UL_i_k*HiB*P_iB_k);
        W_iB_k = inv(E_iB_k_star);
        E_iB_k = eye(fdcomm.ULstream_num(ii))-2*U_iB_k*HiB*P_iB_k+U_iB_k*R_UL_i_k*U_iB_k';
        W_UL{ii,k} = W_iB_k;
        U_UL{ii,k} = U_iB_k;
        E_UL_star{ii,k} = E_iB_k_star;
        E_UL{ii,k} = E_iB_k;
    end
    for jj = 1 : J
        HBj = fdcomm.DLchannels{jj};
        R_DL_j_k = cov.total_DL{jj,k};
        P_Bj_k = fdcomm.DLprecoders{jj,k};
        U_Bj_k = P_Bj_k'*HBj'/(R_DL_j_k);
        E_Bj_k_star = abs(eye(fdcomm.DLstream_num(jj))-P_Bj_k'*HBj'/R_DL_j_k*HBj*P_Bj_k);
        W_Bj_k = inv(E_Bj_k_star);
        E_Bj_k = abs(eye(fdcomm.DLstream_num(jj))-2*U_Bj_k*HBj*P_Bj_k+U_Bj_k*R_DL_j_k*U_Bj_k');
        W_DL{jj,k} = W_Bj_k;
        U_DL{jj,k} = U_Bj_k;
        E_DL_star{jj,k} = E_Bj_k_star;
        E_DL{jj,k} = E_Bj_k;
    end
end
fdcomm.UL_WMMSE_RX = U_UL;
fdcomm.DL_WMMSE_RX = U_DL;
fdcomm.UL_weights = W_UL;
fdcomm.DL_weights = W_DL;
fdcomm.UL_MMSE = E_UL_star;
fdcomm.DL_MMSE = E_DL_star;
fdcomm.UL_MMSE_nop = E_UL;
fdcomm.DL_MMSE_nop = E_DL;
fdcomm.UL_rate = 
end

