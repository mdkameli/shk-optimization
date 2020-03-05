l = 10;
pd=[];ga=0;v=[];
G2 = G(1:4*K,1:4*K);
u = diag(G(1:4*K,1:4*K));
sig = G2 - u*u';
for i=1:l
    pd = [pd normrnd(u,diag(sig))];
    for j=1:4*K
        ga = qfuncinv(u(j))*sqrt(sig(j,j))+u(j);
        v(j,i) = (sign(pd(j,i)-ga)+1)/2;
        ga=0;
    end
end
    
    