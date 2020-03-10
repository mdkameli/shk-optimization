close all; clear all; clc; 
%%%% FOR K Consequtive Task
%% Define Variables
K = 10;                                 % Number of tasks
Tmax = 4;                               % service threshold
d_k = unifrnd(200,800,K,1);
%w_k = 30;                               % CPU cycle for each bit in task (cycl/bit)
eng = []; d_ave = []; FT_ave = []; w_ave = []; opt_eng_tmp = 0;
opt_eng=[]; eng_loc=[]; eng_rel=[]; eng_redg=[]; eng_edg=[];
%% Random Graph Generation (for K nodes regarding K tasks)
A = graph_generator(K);
grap = simplify(digraph(A));
plot(grap)
%% Problem Matrix Generation
for i=1:30
    for j=1:1
        %d_k = unifrnd(50+100*i,150+100*i,K,1);
        w_k = 1.5*i;
        [M0, Mj, Mkp, Mkd, Mkrj, Mkr , b2] = Mat_Gen(d_k, w_k, K);
        % Optimization Formulation
        G = cvx_opt(M0, Mj, Mkp, Mkd, Mkrj, Mkr, grap, A, K, Tmax);
        % If G is not of rank 1
        while 1
            v = v_formulation(G, K);
            % Finding FTk and RTk but in this case we also need the directed acyclic graph which is denoted here as 'grap'
            FT = find_FT(b2, v, grap, Tmax, K);
            if FT(K) > Tmax
                continue
            else
                break
            end
        end
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

for i=1:30
    w_k = 1.5*i;
    [M0, Mj, Mkp, Mkd, Mkrj, Mkr , b2] = Mat_Gen(d_k, w_k, K);
    % LOCAL
    v_loc = different_v(1, K);
    FT_loc = find_FT(b2, v_loc, grap, Tmax, K);
    sol_loc = [v_loc FT_loc.' 1 1];
    energy_loc = sol_loc*M0*sol_loc.';
    eng_loc = [eng_loc energy_loc];
    % Device to Relay
    v_rel = different_v(2, K);
    FT_rel = find_FT(b2, v_rel, grap, Tmax, K);
    sol_rel = [v_rel FT_rel.' 1 1];
    energy_rel = sol_rel*M0*sol_rel.';
    eng_rel = [eng_rel energy_rel];
    % Relay to Edge
    v_redg = different_v(3, K);
    FT_redg = find_FT(b2, v_redg, grap, Tmax, K);
    sol_redg = [v_redg FT_redg.' 1 1];
    energy_redg = sol_redg*M0*sol_redg.';
    eng_redg = [eng_redg energy_redg];
    % Device to Edge
    v_edg = different_v(4, K);
    FT_edg = find_FT(b2, v_edg, grap, Tmax, K);
    sol_edg = [v_edg FT_edg.' 1 1];
    energy_edg = sol_edg*M0*sol_edg.';
    eng_edg = [eng_edg energy_edg];
end
%% PLOT THE RESULTS
plot(w_ave, opt_eng, 'r-o')
grid on
hold on
plot(w_ave, eng_loc, '-*')
plot(w_ave, eng_rel, '-+')
plot(w_ave, eng_redg, '-^')
plot(w_ave, eng_edg, '-.>')
legend({'Optimum Solution', 'Local', 'Relay', 'Relay to Edge', 'Device to Edge'}, 'Location','northwest')
hold off 


