function [radar] = radar_MMSE(radar, cov)
% radar_MMSE returns the MMSE receiver and Weight Matrice
%   Detailed explanation goes here
Mr = radar.TX;
eta_T = radar.channelgain;
for nr = 1 : radar.RX
    eta_t = eta_T(nr);
    St_nr = cov.S_tr(:,:,nr);
    R_in_nr = cov.inr(:,:,nr);
    Urnr= eta_t^2*St_nr'/(eta_t^2*(St_nr*St_nr')+R_in_nr);
    Ernr_star = abs(eta_t^2*(eye(Mr)-eta_t^2*St_nr'/(R_in_nr+St_nr*St_nr')*St_nr));
    radar.WMMSE_weights{nr,1} = inv(Ernr_star);
    Ernr = abs(eta_t^2*(eye(Mr)-2*Urnr*(St_nr*St_nr')*Urnr')+Urnr*R_in_nr*Urnr');
    radar.WMMSE_RX{nr,1}= Urnr;
    radar.MMSE{nr,1} = Ernr_star;
    radar.MMSE_nop{nr,1} = Ernr;
end

