function A = graph_generator(K)
%% Random Graph Generation (for K nodes regarding K tasks)
while 1
    A = round(rand(K));
    A = triu(A) + triu(A,1)';
    A = A - diag(diag(A));
    A = triu(A);
    if length(find(all(A == 0,2))) > 1
        continue
    else
        break
    end
end
end



