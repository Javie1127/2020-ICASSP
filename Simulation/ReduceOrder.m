function [Rt,bt,jeq,jieq,Anew,bnew] = ReduceOrder(Aeq,beq,ub,lb,A,b)
%Reduce number of unknowns via Aeq*[x_eq ; x_ieq] =beq
%check inputs
N =  size(Aeq,2); %number of independent variables
assert(size(Aeq,1) == size(beq,1))%size should match
assert(size(beq,2) == 1) %column vector only
assert(size(Aeq,2)>size(Aeq,1)) % Aeq should be under-determint.
assert(size(ub,1) == N)%shape requirements
assert(size(lb,1) == N)
assert(size(ub,2) == 1)
assert(size(lb,2) == 1)
assert(size(A,1) == size(b,1))%size should match
assert(size(b,2) == 1) %column vector only
assert(size(A,2) == N)

[R,jeq]=rref([Aeq,beq]);%Gauss eliminiation
%x(jeq) = x_eq
%the rest is x_ieq
%x_eq = bt - Rt*x_ieq;
jieq = setdiff(1:N,jeq);%get the index of x_ieq
Rt = R(:,jieq);
bt = R(jeq,end);
%generate new ieq constraints matrix
Anew1 = A(:,jieq)-A(:,jeq)*Rt;
bnew1 = b - A(:,jeq)*bt;
Anew2 = -Rt;
bnew2 = ub(jeq)-bt;
Anew3 = Rt;
bnew3 = bt-lb(jeq);

Anew = [Anew1;Anew2;Anew3;eye(numel(jieq));-eye(numel(jieq))];
bnew = [bnew1;bnew2;bnew3;ub(jieq);-lb(jieq)];
end


