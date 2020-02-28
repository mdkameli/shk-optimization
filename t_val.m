function [upp] = t_val(v)
%function to obtain value for t-distributions with v degrees of freedom for
%two-sided critical regions and 95% confidence interval
positive_side = 1-0.05/2;   
negative_side = 0.05/2;
upp = tinv(positive_side,v);
low = tinv(negative_side,v);
end

