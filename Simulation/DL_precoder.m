function [P_Bj_k] = DL_precoder(jj, k,...
    H_DL, H_BB, W_DL_k, W_UL_k, Wr, U_DL_k, U_UL_k,Ur,...
    eta_B, eta_t, d_DL_n, d_DL_1, P_DL_k, Ar, f_Btr,F_rtr,...,
    lambda_DL,mu_DL,mu_UL, napla_PBj_R_DL_k,napla_PBj_R_UL_k)
% DL_precoder returns $P^star_{i,B}[k]$
%-----------------------------------
% W_DL_k: dj*dj*J DL UE Weight Matices
% W_UL_k: di*di*I UL UE Weight Matrices
% lambda_DL: Lagrange multiplier vector for the DL power constraints
% mu_DL: K*J Lagrange multiplier vector for the DL rate constraints
% mu_UL: K*I Lagrange multiplier vector for the UL rate constraints
% napla_PBj_R_DL_k: Mc*dj*J
% napla_PBj_R_UL_k: Mc*dj*I
%%
Mr = size(H_DL,2);
J = size(H_DL,3);
I = size(W_UL_k,3);
for g = 1:J
 H_Bg = H_DL{g};
 W_Bg_k = W_DL_k{g};
 U_Bg_k = U_DL_k{g};
 A_temp_1 = H_Bg'*U_Bg_k'*W_Bg_k*U_Bg_k*H_Bg;
 A_1 = A_temp_1+A_1;
 mu_gd_k = mu_DL(k,g);
 napla_PBj_RBg = napla_PBj_R_DL_k(:,:,g);
 temp_rate_DL = mu_gd_k*napla_PBj_RBg+temp_rate_DL;
end
for q =1:I
    U_iB_k = U_UL_k{q};
    W_iB_k = W_UL_k{q};
    A_temp_2 = H_BB'*U_iB_k'*W_iB_k*Ui_B_k*H_BB;
    A_2 = A_2+A_temp_2;
    mu_qu_k = mu_UL(k,q);
    napla_PBj_RqB = napla_PBj_R_UL_k{q};
    temp_rate_UL = mu_qu_k*napla_PBj_RqB+temp_rate_UL;
end
lambda_DL_k = lambda_DL(k);
A = A_1+A_2+lambda_DL_k*eye(Mr);
d_Bj_k_n = d_DL_n(:,jj,k);
d_Bj_k_1 = d_DL_1(:,jj,k);
for nr = 1: Nr
    Urnr = Ur(:,:,nr);
    urnr_k = Urnr(:,k);
    Wrnr = Wr(:,:,nr);
    B_temp = eta_B^2*d_Bj_k_n* urnr_k'*Wrnr*urnr_k*d_Bj_k_n'+...
       eta_t^2*d_Bj_k_1*urnr_k'*Wrnr*urnr_k*Wrnr*urnr_k*d_Bj_k_1;
   B = B_temp+B;
end
%% Matrix C
CC = HBj'*UBj_k'*WBj_k+temp_rate_DL*temp_rate_UL;
for nr = 1 : Nr
    Urnr = Ur(:,:,nr);
    urnr_k = Urnr(:,k);
    Wrnr = Wr(:,:,nr);
    f_Btnr = f_Btr(nr);
    kk = 0:1:K-1;
    f_tnr = F_rtr(:,nr);
    q_Btnr = exp(1i*2*pi*kk*f_Btr);
    Q_Btnr = diag(q_Btnr);
    q_Btnr_k = q_Btnr(k);
    q_rnr_k = exp(1i*2*pi*(k-1)*f_tnr);
    Qrnr_k = diag(q_rnr_k);
    C1 = eta_t^2*(q_Btnr_k-cos(2*pi*(k-1)*f_Btnr)*Qrnr_k*ar_k*urnr_k')*Wrnr*urnr_k*d_Bj_k_1+C1;
    for g =1:J
        while g ~= jj
            P_Bg_k = P_DL_k{g};
            d_Bg_k_n = d_DL_n{g};
            C2_temp_1 = P_Bg * d_Bg_k_n*urnr_k'*Wrnr*urnr_k*d_Bj_k_n+C2_temp_1;
            d_Bg_k_1 = d_DL_1(:,g,k);
            C2_temp_2 = P_Bg*d_Bg_k_1*urnr_k'*Wrnr*urnr_k*d_Bg_k_1'+C2_temp_2;
        end
    end
    C2 = eta_B^2*C2_temp_1+eta_t^2*C2_temp_2+C2;
    for m = 1:K
        while m ~= K
            urnr_m = Urnr(L,m);
            a_m = Ar(m,:);
            q_rnr_m = exp(1i*2*pi*(m-1)*f_tnr);
             Qrnr_m = diag(q_rnr_m);
            for g =1:J
                P_Bg_m = P_DL_k(:,:,g,m);
                d_Bg_m_n = d_DL_n(:,g,m);
                d_Bg_m_1 = d_DL_1(:,g,m);
                temp1 = P_Bg_m*d_Bg_m_1+temp1;
                temp2 = P_Bg_m*d_Bg_m_n+temp2;
            end
            C3_temp_1 = (Qrnr_m*a_m+temp1)*urnr_m'*Wrnr*urnr_k*d_Bj_k_1+C3_temp_1;
            C3_temp_2 = temp2*urnr_m'*Wrnr*urnr_k*d_Bj_k_n+C3_temp_2;
        end
    end
    C3 = eta_t^2*C3_temp_1+eta_B^2*C3_temp_2+C3;
end
HBj = H_DL(:,:,jj);
UBj_k = U_DL_k(:,:,jj);
WBj_k = W_DL_k(:,:,jj);
C = CC+C1-C2-C3;
P_Bj_k = sylvester(A,B,C);
end
