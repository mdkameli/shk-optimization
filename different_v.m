function v = different_v(a, K)
%%Function is for assuming v when all the tasks computation occurs at only one of these locations: local, relay, edge via relay and direct edge

% Optimal Solution
disp('Optimal Solution');      % For Computed v;
% 1) All computation done locally
if a == 1    
    v = repmat([1 0 0 0], 1, K);
    disp('All computation done locally');
% 2) All computation done in Relay
elseif a == 2    
    v = repmat([0 1 0 0], 1, K);
    disp('All computation done in Relay');
% 3) All computation done in edge through relay connection
elseif a == 3  
    v = repmat([0 0 1 0], 1, K);
    disp('All computation done in edge through relay connection');
% 4) All computation done in edge through device direct connection
else
    v = repmat([0 0 0 1], 1, K);
    disp('All computation done in edge through device direct connection');
end
end

