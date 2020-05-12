close all; clear all; clc; 
%%%% FOR K Consequtive Task
%% Define Variables
K = 15;                                 % Number of tasks
Tmax = 3;                               % service threshold
rep = 100;                               % Solution repetition
rang = 16;                              % Solution range of change
d_k = unifrnd(300, 500, K, 1);
%w_k = 30;                               % CPU cycle for each bit in task (cycl/bit)
eng = []; d_ave = []; FT_ave = []; w_ave = [];% opt_eng_tmp = 0;
eng_loc=[]; eng_rel=[]; eng_redg=[]; eng_edg=[];
eng_loc_j=[]; eng_rel_j=[]; eng_redg_j=[]; eng_edg_j=[];
%% Random Graph Generation (for K nodes regarding K tasks)
A = graph_generator(K);
grap = simplify(digraph(A));
figure 
plot(grap)
%% Problem Matrix Generation
for i=1:rang
    for j=1:rep
        %d_k = unifrnd(50+100*i,150+100*i,K,1);
        w_k = 5*i;
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
    %opt_eng_tmp = mean(eng);
    %eng=[];
    %opt_eng = [opt_eng opt_eng_tmp];
    %fprintf ("Energy Consumption = %f (j)\n", energy_consumption)
end
%% Different Scenarios: along with K input 1 for local, 2 for relay, 3 for edge via relay and 4 for edge execution in following function
%v = different_v(1,K); %here as example input as 1 will give local execution only
for i=1:16
    w_k = 5*i;
    [M0, Mj, Mkp, Mkd, Mkrj, Mkr , b2] = Mat_Gen(d_k, w_k, K);
    % LOCAL
    v_loc = different_v(1, K);
    FT_loc = find_FT(b2, v_loc, grap, Tmax, K);
    sol_loc = [v_loc FT_loc.' 1 1];
    energy_loc = sol_loc*M0*sol_loc.';
    eng_loc = [eng_loc energy_loc];
    eng_loc_j = [eng_loc_j energy_loc*FT_loc(K)];
    % Device to Relay
    v_rel = different_v(2, K);
    FT_rel = find_FT(b2, v_rel, grap, Tmax, K);
    sol_rel = [v_rel FT_rel.' 1 1];
    energy_rel = sol_rel*M0*sol_rel.';
    eng_rel = [eng_rel energy_rel];
    eng_rel_j = [eng_rel_j energy_rel*FT_rel(K)];
    % Relay to Edge
    v_redg = different_v(3, K);
    FT_redg = find_FT(b2, v_redg, grap, Tmax, K);
    sol_redg = [v_redg FT_redg.' 1 1];
    energy_redg = sol_redg*M0*sol_redg.';
    eng_redg = [eng_redg energy_redg];
    eng_redg_j = [eng_redg_j energy_redg*FT_redg(K)];
    % Device to Edge
    v_edg = different_v(4, K);
    FT_edg = find_FT(b2, v_edg, grap, Tmax, K);
    sol_edg = [v_edg FT_edg.' 1 1];
    energy_edg = sol_edg*M0*sol_edg.';
    eng_edg = [eng_edg energy_edg];
    eng_edg_j = [eng_edg_j energy_edg*FT_edg(K)];
end
ax=[eng_loc;eng_rel;eng_redg;eng_edg];
ay=min(ax);
%% Average Time & Energy of Optimum soloution
FT_vec =[]; opt_eng=[];
for i=1:rep:rep*rang
    FT_vec = [FT_vec mean(FT_ave(i:i+rep-1))];
    opt_eng = [opt_eng mean(eng(i:i+rep-1))];
end
opt_eng_j = opt_eng.*FT_vec;
%% PLOT THE RESULTS IN WATT
figure 
plot(w_ave, opt_eng, 'r-d','MarkerSize',8,'LineWidth',1.5)
grid on
hold on
plot(w_ave, eng_loc, '--*')
plot(w_ave, eng_rel, '--+')
plot(w_ave, eng_redg, '--^')
plot(w_ave, eng_edg, '-->')
%plot(w_ave, ay, 'b-s','MarkerSize',10,'LineWidth',1.5)
legend({'Optimum Solution', 'Local', 'Relay', 'Relay to Edge', 'Device to Edge', 'Base Solution'}, 'Location','northwest')
hold off 
%% PLOT THE RESULTS IN JULE
figure 
plot(w_ave, opt_eng_j, 'r-d','MarkerSize',8,'LineWidth',1.5)
grid on
hold on
plot(w_ave, eng_loc_j, '--*')
plot(w_ave, eng_rel_j, '--+')
plot(w_ave, eng_redg_j, '--^')
plot(w_ave, eng_edg_j, '-->')
%plot(w_ave, ay, 'b-s','MarkerSize',10,'LineWidth',1.5)
legend({'Optimum Solution', 'Local', 'Relay', 'Relay to Edge', 'Device to Edge'}, 'Location','northwest')
hold off 
%% RESET
%eng_loc=[];eng_rel=[];eng_redg=[];eng_edg=[];
%eng_loc_j=[];eng_rel_j=[];eng_redg_j=[];eng_edg_j=[];









