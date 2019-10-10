%%% Convergence analysis

[fdcomm,radar,radar_comm] = parameters_icassp;
%% Initialization
%[fdcomm, radar] = ini_icassp(fdcomm, radar);
%% WMMSE algorithm
[fdcomm, radar, radar_comm,Xi_mse] = WMMSE_algorithm_ICASSP(fdcomm, radar, radar_comm);
%% Plot the convergence 
iter_max = length(Xi_mse);
iter = 1:iter_max;