function [a_mr_bigstar] = Nearest_PAR(radar)
%function [a_mr_bigstar] = Nearest_PAR(z,rho,c)
%Nearest_PAR returns the nearest a 
%
%% Initialization
rho = radar.rho;
d = size(z,1);
z_nor = normalize(z,'norm',2);
sigma = sqrt(c*rho/d);
k = 0;
a_mr_bigstar = zeros(d,1);
a_mr_star_mag = abs(z_nor); % finding the magnitude
%% 
M = find(a_mr_star_mag == min(a_mr_star_mag));
k = d - size(M,1);
%%%%%%%%%% Judge if the set is uniquely defined
%while size(z_nor(M),1) > size(unique(z_nor(M)),1)
%    k=k+1;
%   M = find(a_mr_star_mag(1:d-k) == min(a_mr_star_mag(1:d-k)));
%end
z_min = min(a_mr_star_mag);
if z_min == 0
    for i = 1:d
        if ismember(i,M) 
            a_mr_bigstar(i) = sqrt((c-k*sigma^2)/(d-k));
        else
            a_mr_bigstar(i) = sigma*exp(1i*angle(z_nor(i)));
        end
    end
else
    gamma=sqrt((c-k*sigma^2)/norm(z_nor(M))^2);
    while gamma*z_nor(M)>sigma
        k=k+1;
        M = find(a_mr_star_mag(1:d-k) == min(a_mr_star_mag(1:d-k)));
    end
    for i = 1:d
        if ismember(i,M) 
            a_mr_bigstar(i) = gamma*z_nor(i);
        else
            a_mr_bigstar(i) = sigma*exp(1i*angle(z_nor(i)));
        end
    end
end
end

