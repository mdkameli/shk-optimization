function v = v_formulation(G, K)
%When rank of G matrix is not 1 then find vector v which is the offloading
%decision
%rng(3);
e3=@(j) [zeros(j-1,1);1;zeros(4-j,1)].';
zeta = unifrnd(0,1,K,1);
gama = G(:,end);
gama = gama(1:end-K-2);
v = [];
for j=1:K
    if zeta(j) <= gama((j-1)*4+1)
        v = [v e3(1)];
    elseif gama((j-1)*4+1) < zeta(j) <= gama((j-1)*4+2)
        v = [v e3(2)];
    elseif gama((j-1)*4+2) < zeta(j) <= gama((j-1)*4+3)
        v = [v e3(3)];
    else
        v = [v e3(4)];
    end
end
end

