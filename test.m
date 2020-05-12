clear all; close all;
% Define Variables
K = 15;                                  % Number of tasks
Tmax = 3;  
d_k = unifrnd(300,500,K,1);             % data size (uniformly distributed btw 300-500 KB)
w_k = 10;                               % CPU cycle for each bit in task (cycl/bit)
floc_k = 0.5;                           % local CPU cycl frequency (uniform dist. btw 0.1-0.5 Gcyc/s)
fedg_k = 2;                             % edge CPU cycle frequency (Gcyc/s)
frel_k = 1;                             % assume 2nd device(relay) CPU cycle frequency (Gcyc/s)
k1 = 1e-27;                             % Effective switched capacitance
ptra_k = 2*1e-3;                           % Idle transmission power (mW) (of local node)
pcir_k = repmat(0.01, K, 1);             % Idle circuit power (uniform dist. btw 0.001-0.01 w)
ptra_rel_k = 4*1e-3;                       % Transmission power (mW) of relay
%%%% CHANNEL PARAMETER DEFINITION
% CHANELS NOISES
sig = 1e-9;                         % Varriance of the Gaussian channel Noise
% CHANNELS GAINS
H_dev_k = 1e-7;                         % Signal Energy (channel gain between the node and relay(2nd device))
H_rel_k = 1e-6;                     % channel gain between the relay and edge
H_edg_k = 1e-7;                     % channel gain between the node and edge
% CHANNELS BANDWIDTH
B = 5; %randi([3,6], K, 1);             % Channel Bandwidth (Mhz)
%% Initial Vector Computation
d_k = 8*1000.*d_k;                                                                          %% KB to bit
C_k = w_k.*d_k;                                                                             % CPU required to complete task k (cycle)
%%%% CHANNELS RATES
B = B*1e6;                                                                          %% Mhz to hz
r_k = B*log2(1+(ptra_k*H_dev_k/sig));                                               % Transmission rate btw device and relay (bit/s)
r_rel_k = B*log2(1+(ptra_rel_k*H_rel_k/sig));                                       % Transmission rate btw relay and edge (bit/s)
r_edg_k = B*log2(1+(ptra_k*H_edg_k/sig));                                           % Transmission rate btw device and edge (bit/s)
%%%% TIME & ENERGY COMPUTATION
floc_k = 1e9.*floc_k;                                                                       %% Gcyc/s to cyc/s
frel_k = 1e9.*frel_k;                                                                       %% Gcyc/s to cyc/s
fedg_k = 1e9.*fedg_k;                                                                       %% Gcyc/s to cyc/s
% LOCAL %
Tloc_k = C_k./floc_k;                                                                       % Local transmission time
Eloc_k = (k1*floc_k.^2).*C_k;                                                               % Local energy cunsumption
% DEVICE TO RELAY
Td2d_k = d_k./r_k + C_k./frel_k;                                                            % Transmission time when data send from local node --> relay
Ed2d_k = ptra_k*(d_k./r_k) + pcir_k.*(C_k./frel_k);                                         % Energy consumption when execution done at relay
% RELAY TO EDGE %
Trel_edg_k = d_k./r_k + d_k./r_rel_k + C_k./fedg_k;                                         % Transmission time when data send from local node --> relay --> edge
Erel_edg_k = ptra_k*(d_k./r_k) + pcir_k.*(d_k./r_rel_k + C_k./fedg_k);                      % Energy consumption when data send from local node  --> relay --> edge
% DEVICE TO EDGE %
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