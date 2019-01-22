function [N,k1,k2,T1,T2] = find_darboux_frame(P,PS)
% An algorithm to estimate the darboux frame (normal, principal directions and principal 
% curvatures) when given a point and its neighbours.
% This algorithm was developed by Eyal Hameiri and Ilan Shimshoni
% For details please refer to 3DPVT 2002.
for i=1:size(PS,1)
    ps2(i,1:3) = PS(i,1:3)-P;
end;

[U,S,V] = svd(ps2);
N = V(:,3);
if(N(3) < 0)
    N= -N;
end;
Mv = zeros(3,3);
A=0;
B=0;
C=0;
for i=1:size(PS,1)
    Kp(i) = 2*N'*ps2(i,:)'/norm(ps2(i,:))^2;
    T(i,1:3) = ps2(i,:)-ps2(i,:)*N*N';
    T(i,:) = T(i,:)/norm(T(i,:))';
    theta(i) = asin(max(min(T(i,:)*T(1,:)',1),-1));
    Mv = Mv + Kp(i) * T(i,:)'*T(i,:);
    A=A+cos(theta(i))^4;
    B=B+(cos(theta(i))*sin(theta(i)))^2;
   C=C+sin(theta(i))^4;
end;
[U,S,V]=svd(Mv);
e1=S(1,1);
e2=S(2,2);
m= [A,B
    B,C];
k12=inv(m)*[e1,e2]';
if(k12(1) > k12(2))
    k1 = k12(1);
    k2 = k12(2);
    T1 = V(:,1);
    T2 = V(:,2);
else
    k1 = k12(2);
    k2 = k12(1);
    T1 = V(:,2);
    T2 = V(:,1);
end;
