close all; clear all; clc; 
%%%% FOR K Consequtive Task
%% Define Variables
rng(1);
K = 15;                                 % Number of tasks
Tmax = 10;                               % service threshold
d_k = unifrnd(300,500,K,1);
%w_k = 30;                               % CPU cycle for each bit in task (cycl/bit)
eng = []; d_ave = []; FT_ave = []; w_ave = [];opt_eng_tmp = 0; opt_eng=[];
%% Random Graph Generation (for K nodes regarding K tasks)
A = graph_generator(K);
grap = simplify(digraph(A));
plot(grap)
%% Problem Matrix Generation
for i=1:5                   %5 different values of w_k taken and find its respective optimized energy
    %d_k = unifrnd(50+100*i,150+100*i,K,1);
    w_k = 5*i;
    [M0, Mj, Mkp, Mkd, Mkrj, Mkr , b2] = Mat_Gen(d_k, w_k, K);
    for j=1:5             %L = 5 times the samples are run to find optimized energy
        %% Optimization Formulation
        G = cvx_opt(M0, Mj, Mkp, Mkd, Mkrj, Mkr, grap, A, K, Tmax);
        %% If G is not of rank 1
        v = v_formulation(G, K);
        %% Finding FTk and RTk but in this case we also need the directed acyclic graph which is denoted here as 'grap'
        FT = find_FT(b2, v, grap, Tmax, K);
        % Final Solution
        opt_sol = [v FT.' 1 1];
        energy_consumption = opt_sol*M0*opt_sol.';
        eng = [eng energy_consumption];
        %d_ave = [d_ave mean(d_k)];
        FT_ave = [FT_ave FT(K)];        
        %fprintf ("Energy Consumption = %f (j)\n", energy_consumption)
        %fprintf ("Average Data Size  = %f (KB) ,i=%d\n", mean(d_k), i)        
    end
    w_ave = [w_ave w_k];
    opt_eng_tmp = min(eng);
    eng=[];
    opt_eng = [opt_eng opt_eng_tmp];
    %fprintf ("Energy Consumption = %f (j)\n", energy_consumption)
end
%% Different Scenarios: along with K input 1 for local, 2 for relay, 3 for edge via relay and 4 for edge execution in following function
%v = different_v(1,K); %here as example input as 1 will give local execution only
plot(w_ave, opt_eng, 'r-o')
grid on
