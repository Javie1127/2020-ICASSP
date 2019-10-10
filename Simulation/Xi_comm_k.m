function [fdcomm] = Xi_comm_k(fdcomm, k)
Xi_UL_k = 0;
I = fdcomm.UL_num;
J = fdcomm.DL_num;
for ii = 1:I
    alpha_iu = fdcomm.alpha_UL(ii);
    W_iu_k = fdcomm.UL_weights{ii,k};
    E_iu_k = fdcomm.UL_MMSE_nop{ii,k};
    Xi_UL_k = Xi_UL_k + alpha_iu*abs(trace(W_iu_k*E_iu_k));
end
Xi_DL_k = 0;
for jj = 1:J
    alpha_jd = fdcomm.alpha_DL(jj);
    W_jd_k = fdcomm.DL_weights{jj,k};
    E_jd_k = fdcomm.DL_MMSE_nop{jj,k};
    Xi_DL_k = Xi_DL_k+alpha_jd*abs(trace(W_jd_k*E_jd_k));
end
fdcomm.Xi_UL(k) = Xi_UL_k;
fdcomm.Xi_DL(k) = Xi_DL_k;
end

