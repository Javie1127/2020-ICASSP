function [fdcomm, radar] = altermating_projection...
(fdcomm, radar, radar_comm, PAR, iter1_max)

%% Initialization
iter1                   = 0; % iteration index
while iter1< iter1_max
    % P_{i,\B}, P_{\B,j}
    [fdcomm, radar] = WMMSE_algorithm(fdcomm, radar, radar_comm);
    % Calculate A^\star
    radar = Nearest_PAR(radar, PAR);
    iter1 = iter1+1;
end
% Caculate the optimal WMMSE RXs
fdcomm = Comm_MMSE(fdcomm, cov);
radar = radar_MMSE(radar, cov);
end

