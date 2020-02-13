close all; clear all; clc; 
%%%% FOR K Consequtive Task
%% Define Variables
K = 10;                                 % Number of tasks
B = 5;                                  % Bandwidth (Mhz)
sig = 1e-9;                             % Noise Energy (w)
H_k = 1e-6;                             % Signal Energy
d_k = unifrnd(300,500,K,1);             % data size (uniformly distributed btw 300-500 KB)
w_k = 30;                               % CPU cycl for each bit in task (cycl/bit)
Tmax = 4;                               % service threshold
floc_k = unifrnd(0.1,0.5,K,1);          % local CPU cycl frequency (uniform dist. btw 0.1-0.5 Gcyc/s)
fedg_k = 2;                             % edge CPU cycl frequency (Gcyc/s)
fclo_k = 4;                             % Cloud CPU cycl frequency (Gcyc/s)
k1 = 1e-27;                             % Effective switched capacitance
ptra_k = 0.1;                           % Idle transmission power (w)
pcir_k = unifrnd(0.001,0.01,K,1);       % Idle circuit power (uniform dist. btw 0.001-0.01 w)
R_EC = 5;                               % Transmission rate btw relay and cloud (MB/s)
%% Initial Vector Computation
d_k = 8*1000.*d_k;                                              %% KB to bit
C_k = w_k.*d_k;                                                 % CPU required to complete task k (cycle)
B = B*1e6;                                                      %% Mhz to hz
r_k = B*log2(1+(ptra_k*H_k/sig));                               % Transmission rate btw device and relay (bit/s)
floc_k = 1e9.*floc_k;                                           %% Gcyc/s to cyc/s
Tloc_k = C_k./floc_k;                                           % Local transmission time
Eloc_k = (k1*floc_k.^2).*C_k;                                    % Local energy cunsumption
fedg_k = 1e9*fedg_k;                                            %% Gcyc/s to cyc/s
Tedg_k = d_k./r_k + C_k./fedg_k;                                % Relay transmission time
Eedg_k = ptra_k*(d_k./r_k) + pcir_k.*(C_k./fedg_k);             % Realy energy consumption
fclo_k = 1e9.*fclo_k;                                           %% Gcyc/s to cyc/s
R_EC = 1e6*8*R_EC;                                              %% MB/s to bit/s
Tclo_k = d_k./r_k + d_k./R_EC + C_k./fclo_k;                    % Cloud transmission time
Eclo_k = ptra_k*(d_k./r_k) + pcir_k.*(d_k./R_EC + C_k./fclo_k); % Cloud energy consumption
%% Initial Vector re-ordering
b0 = [Eloc_k.';Eedg_k.';Eclo_k.'];
b0 = b0(:);
b0 = [b0 ; pcir_k ; 0];
b1= [zeros(1,3*K-3), Tloc_k(K), Tedg_k(K), Tclo_k(K), zeros(1,K-1), 1, 0].';
b2 = [Tloc_k.' ;Tedg_k.' ; Eclo_k.'];
b2 = b2(:);
b2 = [b2 ; ones(K,1)];
%% 
e=@(j) [zeros(j-1,1);1;zeros(4*K+1-j,1)];
e2=@(j) [zeros(j-1,1);1;zeros(4*K-j,1)];
bkp=@(k) (e(3*k-2)+e(3*k-1)+e(3*k));
bj=@(j) (e2(3*j-2)+e2(3*j-1)+e2(3*j)+e2(3*K+j));

%%% handel definition
a0=zeros(4*K+1,4*K+1);
a00=zeros(4*K,4*K);
%a1=@(j) (-1/2)*(b2.'*diag(bj(j))).';
%a2=@(k) (1/2)*e(3*K+k);

%% Problem matrix definition
M0=[a0, (1/2)*b0; (1/2)*b0.', 0];
Mj=@(j) [diag(e(j)), (-1/2)*e(j); (-1/2)*e(j).', 0];
Mkp=@(k) [a0, (1/2)*bkp(k); (1/2)*bkp(k).', 0];
Mkd=[a0, (1/2)*b1; (1/2)*b1.', 0];
Mkr=@(k) [a0, (1/2)*e(3*K+k); (1/2)*e(3*K+k).', 0];
Mkrj=@(j,k) [a00, (-1/2)*(b2.'*diag(bj(j))).', (1/2)*e2(3*K+k);...
    (-1/2)*(b2.'*diag(bj(j))), 0, 0; (1/2)*e2(3*K+k).', 0, 0];
%% Optimization Formulation
cvx_solver sedumi
cvx_begin sdp
    variable G(4*K+2, 4*K+2) semidefinite
    minimize (trace(M0*G))
    subject to
%          G == semidefinite(4*K+2)
        for i=1:3*K
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
        trace(Mkr(1)*G) == 0
%%%%%%%%
        G(4*K+1,4*K+1) == 1
        G(4*K+1,4*K+2) == 1
        G(4*K+2,4*K+1) == 1
        G(4*K+2,4*K+2) == 1
cvx_end
G = full(G);
x=sqrt(diag(G));
plot(x(1:3*K),'*');
%% Result
fprintf("Minimum Energy = %f",trace(M0*G));





