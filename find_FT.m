function FT = find_FT(b2, v, grap, Tmax, K)
%Find finish time by calculating ready time of the predecessor tasks and
%check its feasibility
RT = zeros(K,1);    %initalize the Ready time vector as 0s
FT = zeros(K,1);    %initalize the Finish time vector as 0s
T = b2(1:4*K)'.*v;  %Obtain computation time vector as 0s which is currently size 4*K:1
T = T(T~=0);        %To make the size K:1 and only obtain the actual execution time computed at either local, D2D, rel_edg or edg
for i=1:K
     preIDs = predecessors(grap,i);
     if isempty(preIDs)
         RT(i) = 0;
         FT(i) = T(i);
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
    fprintf ("Finish Time        = %f (Sec)\n", FT(K))
else
    disp('not feasible solution')
    fprintf ("Finish Time        = %f (Sec)\n", FT(K))
end
end


