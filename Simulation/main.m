
% Main script of the Algorithm 3
%% Load the parameters
[fdcomm,radar,radar_comm] = parameters;
%% Covariance matrices
cov = covmat(fdcomm, radar, radar_comm);
%% Precodeing matrices and radar code matrix
t = 0;
A = A_ini;
max_iter_1 = 15;
iter_1 = 1;
%% ALternating process begins 
while iter_1 <= max_iter_1
    
    %% radar MMSE receiver 
    
    %% UL & DL MMSE receivers
    
    
    [UiB,UBj,WiB,WBj] = Comm_MMSE(PiB,PBj,HiB,Ru_i,Rd_j);
    [P_UL_star,P_DL_star,A_star] = WMMSE_algorithm(U_UL, U_DL, Ur, W_UL, W_DL,...
            Wr,H_UL,H_DL,HrB,H_r_DL,H_UL_DL,max_iter_1);
    iter_1 = iter_1+1;
end