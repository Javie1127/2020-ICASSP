Aeq = [1,1,1];
Beq = 1;
lb = [0,0,0]';
ub = [1,1,1]';
A = [0,0,1];
b = 0.7;
Npts = 700;
x0 = [0.5,0.2,0.3]';
pts = RandPtsInLinearConstraints(Npts,x0,Aeq,beq,ub,lb,A,b);
scatter3(pts(1,:),pts(2,:),pts(3,:))
grid on
xlim([0,1])
ylim([0,1])
zlim([0,1])