function G = cvx_opt(M0, Mj, Mkp, Mkd, Mkrj, Mkr, grap, A, K, Tmax)
%% Function for optimization formulation and to obtain G matrix which is defines as gg'
cvx_solver sedumi
cvx_begin sdp
    variable G(5*K+2, 5*K+2) semidefinite
    minimize (trace(M0*G))
    subject to
%          G == semidefinite(5*K+2)
        for i=1:4*K
            trace(Mj(i)*G) == 0;
        end
%%%%%%%%        
        for i=1:K
            trace(Mkp(i)*G) == 1;
        end
%%%%%%%%        
        trace(Mkd*G) <= Tmax;
%%%%%%%% Task Dependency (K consequtive task)
        for i=1:K
            preIDs = predecessors(grap,i);
            for j=1:length(preIDs)
                trace(Mkrj(preIDs(j),i)) >= 0;
            end
        end
%%%%%%%% Subjected regarding node-1 as a start point
idx = find(all(A == 0,1));
        for i=1:length(idx)
            trace(Mkr(idx(i))*G) == 0;
        end
%%%%%%%%
        G(5*K+1,5*K+1) == 1;
        G(5*K+1,5*K+2) == 1;
        G(5*K+2,5*K+1) == 1;
        G(5*K+2,5*K+2) == 1;
cvx_end
end

