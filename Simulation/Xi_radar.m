function [radar] = Xi_radar(radar)
Nr = radar.RX;
alpha_r = radar.alpha_r; 
Xi_r = zeros(Nr,1);
for nr = 1:Nr
    W_rnr = radar.WMMSE_weights{nr,1};
    E_rnr = radar.MMSE_nop{nr,1};
    Xi_r(nr) = alpha_r(nr)*abs(trace(W_rnr*E_rnr));
end
radar.Xi_r = Xi_r;
end

