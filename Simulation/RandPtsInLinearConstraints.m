function pts = RandPtsInLinearConstraints(Npts,x0,Aeq,beq,ub,lb,A,b)
%generate Npts points in size(x0,1)-dimension space under constraints:
% Aeq*x = beq
% A*x<=b
% lb<=x<=ub
% input: Npts: number of points to generate
%        x0: initial seed point, column vector
%        Aeq,beq,ub,lb,A,b: linear constraint matrix
% output: double array, size: size(x0,1)*Npts
% alogrithm: random walk, choose a direction randomly, then walk random
% length forward or backward.

%check x0
assert(size(x0,2) == 1);
assert(all(abs(Aeq*x0-beq)<= 100*eps));
assert(all(ub>=x0));
assert(all(lb<=x0));
assert(all(A*x0<=b));

% reduce the number of unknowns
Nx = size(x0,1);
[Rt,bt,jeq,jieq,Anew,bnew] = ReduceOrder(Aeq,beq,ub,lb,A,b);
xr0 = x0(jieq); %x0_ieq
Nr = numel(jieq);% number of x_ieq

%initialization
pts = nan(Nx,Npts);
pointx = nan(Nx,1);

for i=1:Npts    
%generate directions randomly
D = normrnd(0,1,Nr,1);
D = D./sqrt(sum(D.^2));
%solve Anew*(r*D+xr0)<=bnew
%get AD.*xr <= bAx
AD = Anew*D;
bAx =(bnew-Anew*xr0);
%get the limits of the step size (ra, rb)
ra=NaN;
rb=NaN;
for j=1:size(AD,1)
    if(AD(j)==0)% parallel with some constraint
        continue
    elseif(AD(j)>0)
       rb = min(rb,bAx(j)/AD(j)); % min() ignores NaN :-)
       continue
    elseif(AD(j)<0)
       ra = max(ra,bAx(j)/AD(j));
    end
end
% check (ra,rb)
if (ra>rb) 
    error('The constraints are contradictory');
end
if (ra>0 || rb<0) 
    warning('Point jumps out of the feasible region');
end

%generate step size randomly
r = ra+(rb-ra)*rand();

xr0 = (xr0+r*D);
pointx(jieq)=xr0;
pointx(jeq) = bt-Rt*xr0;
pts(:,i) = pointx;
end

end