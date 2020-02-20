close all; clear all; clc; 
%%%% FOR K Consequtive Task
%% Define Variables
K = 10;                                 % Number of tasks
B = 5;                                  % Bandwidth (Mhz)
sig = 1e-9;                             % Noise Energy (w) for transmission btw node and relay
sig_rel = 0.75*1e-9;                    % Noise Energy (w) for transmission btw relay and edge  
sig_edg = 1.2*1e-9;                     % Noise Energy (w) for transmission btw node and edge 
H_k = 1e-6;                             % Signal Energy (channel gain between the node and relay(2nd device))
d_k = unifrnd(300,500,K,1);             % data size (uniformly distributed btw 300-500 KB)
w_k = 30;                               % CPU cycle for each bit in task (cycl/bit)
Tmax = 4;                               % service threshold
floc_k = unifrnd(0.1,0.5,K,1);          % local CPU cycl frequency (uniform dist. btw 0.1-0.5 Gcyc/s)
fedg_k = 2;                             % edge CPU cycle frequency (Gcyc/s)
frel_k = 1;                             % assume 2nd device(relay) CPU cycle frequency (Gcyc/s)
k1 = 1e-27;                             % Effective switched capacitance
ptra_k = 0.1;                           % Idle transmission power (w) (of local node)
pcir_k = unifrnd(0.001,0.01,K,1);       % Idle circuit power (uniform dist. btw 0.001-0.01 w)
ptra_rel_k = 0.2;                       % Transmission power (w) of relay
H_rel_k = 1.5*1e-6;                     % channel gain between the relay and edge
H_edg_k = 1.1*1e-6;                     % channel gain between the node and edge
% B_rel = 8;                            % Relay bandwidth

%% Initial Vector Computation
d_k = 8*1000.*d_k;                                                                          %% KB to bit
C_k = w_k.*d_k;                                                                             % CPU required to complete task k (cycle)
B = B*1e6;                                                                                  %% Mhz to hz
% B_rel = B_rel*1e6;                                                                        %% Mhz to hz
r_k = B*log2(1+(ptra_k*H_k/sig));                                                           % Transmission rate btw device and relay (bit/s)
r_rel_k = B*log2(1+(ptra_rel_k*H_rel_k/sig_rel));                                           % Transmission rate btw relay and edge (bit/s)
r_edg_k = B*log2(1+(ptra_k*H_edg_k/sig_edg));                                               % Transmission rate btw device and edge (bit/s)
floc_k = 1e9.*floc_k;                                                                       %% Gcyc/s to cyc/s
Tloc_k = C_k./floc_k;                                                                       % Local transmission time
Eloc_k = (k1*floc_k.^2).*C_k;                                                               % Local energy cunsumption
frel_k = 1e9.*frel_k;                                                                       %% Gcyc/s to cyc/s
Td2d_k = d_k./r_k + C_k./frel_k;                                                            % Transmission time when data send from local node --> relay
Ed2d_k = ptra_k*(d_k./r_k) + pcir_k.*(C_k./frel_k);                                         % Energy consumption when execution done at relay
fedg_k = 1e9*fedg_k;                                                                        %% Gcyc/s to cyc/s
Trel_edg_k = d_k./r_k + d_k./r_rel_k + C_k./fedg_k;                                         % Transmission time when data send from local node --> relay --> edge
Erel_edg_k = ptra_k*(d_k./r_k) + pcir_k.*(d_k./r_rel_k + C_k./fedg_k);                      % Energy consumption when data send from local node  --> relay --> edge
Tedg_k = d_k./r_edg_k + C_k./fedg_k;                                                        % Transmission time when data send directly from local node --> edge
Eedg_k = ptra_k*(d_k./r_edg_k) + pcir_k.*(C_k./fedg_k);                                     % Energy consumption when data send directly from local node --> edge

%% Initial Vector re-ordering
b0 = [Eloc_k.';Ed2d_k.'; Erel_edg_k.'; Eedg_k.'];
b0 = b0(:);
b0 = [b0 ; pcir_k ; 0];
b1= [zeros(1,4*K-4), Tloc_k(K), Td2d_k(K), Trel_edg_k(K), Tedg_k(K), zeros(1,K-1), 1, 0].';
b2 = [Tloc_k.' ;Td2d_k.' ;Trel_edg_k.' ;Tedg_k.'];
b2 = b2(:);
b2 = [b2 ; ones(K,1)];
%% 
e=@(j) [zeros(j-1,1);1;zeros(5*K+1-j,1)];
e2=@(j) [zeros(j-1,1);1;zeros(5*K-j,1)];
e3=@(j) [zeros(j-1,1);1;zeros(4-j,1)].';
bkp=@(k) (e(4*k-3)+e(4*k-2)+e(4*k-1)+e(4*k));
bj=@(j) (e2(4*j-3)+e2(4*j-2)+e2(4*j-1)+e2(4*j)+e2(4*K+j));

%%% handel definition
a0=zeros(5*K+1,5*K+1);
a00=zeros(5*K,5*K);
%a1=@(j) (-1/2)*(b2.'*diag(bj(j))).';
%a2=@(k) (1/2)*e(4*K+k);

%% Problem matrix definition
M0=[a0, (1/2)*b0; (1/2)*b0.', 0];
Mj=@(j) [diag(e(j)), (-1/2)*e(j); (-1/2)*e(j).', 0];
Mkp=@(k) [a0, (1/2)*bkp(k); (1/2)*bkp(k).', 0];
Mkd=[a0, (1/2)*b1; (1/2)*b1.', 0];
Mkr=@(k) [a0, (1/2)*e(4*K+k); (1/2)*e(4*K+k).', 0];
Mkrj=@(j,k) [a00, (-1/2)*(b2.'*diag(bj(j))).', (1/2)*e2(4*K+k);...
    (-1/2)*(b2.'*diag(bj(j))), 0, 0; (1/2)*e2(4*K+k).', 0, 0];
%% Optimization Formulation
cvx_solver sedumi
cvx_begin sdp
    variable G(5*K+2, 5*K+2) semidefinite
    minimize (trace(M0*G))
    subject to
%          G == semidefinite(5*K+2)
        for i=1:4*K
            trace(Mj(i)*G) == 0
        end
%%%%%%%%        
        for i=1:K
            trace(Mkp(i)*G) == 1
        end
%%%%%%%%        
        trace(Mkd*G) <= Tmax
%%%%%%%% Task Dependency (K consequtive task)
        for i=1:K-1
            trace(Mkrj(i,i+1)) >=0
        end
%%%%%%%% Subjected regarding node-1 as a start point
%            for i=1:K
%                trace(Mkr(i)*G) == 0
%            end
        trace(Mkr(1)*G) == 0
%%%%%%%%
        G(5*K+1,5*K+1) == 1
        G(5*K+1,5*K+2) == 1
        G(5*K+2,5*K+1) == 1
        G(5*K+2,5*K+2) == 1
cvx_end
%G = full(G);
%x=sqrt(diag(G));
%plot(x(1:4*K),'*');
%% Result
%fprintf("Minimum Energy = %f",trace(M0*G));
%% If G is not of rank 1
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
    elseif gama((j-1)*4+3) < zeta(j)
        v = [v e3(4)];
    end
end

%%Finding FTk and RTk but in this case we also need the directed acyclic graph which is denoted here as 'grap'
RT = zeros(K,1);    %initalize the Ready time vector as 0s
FT = zeros(K,1);    %initalize the Finish time vector as 0s
T = b2(1:4*K)'.*v;  %Obtain computation time vector as 0s which is currently size 4*K:1
T = T(T~=0);        %To make the size K:1 and only obtain the actual execution time computed at either local, D2D, rel_edg or edg
for i=1:K
     preIDs = predecessors(grap,i);
     if isempty(preIDs)
         RT(i) = 0;
         FT(i) = T(i)
     else
         tot = numel(preIDs);                      
         max = FT(preIDs(1));
         for j=2:tot
           if FT(preIDs(j)) > max
              max = FT(preIDs(j));  %obtain the highest finish time amongst the predecessor nodes  
           end
         end
             RT(i) = max;
             FT(i) = RT(i) + T(i);
     end
end
if FT(K) <= Tmax
    disp('feasible solution for the optimization problem')
    FT(K)
else
    disp('not feasible solution')
    FT(K)
end
