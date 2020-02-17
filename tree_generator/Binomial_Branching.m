function [ Binomial ] = Binomial_Branching(n,nd,t1,p0)
%Zero branching is allowed after the second generation thus the branching at the root is (1-nd). 
%For the rest, it is (0-nd). The branching is generated using Binomial random numbers with p=p0, and n=nd, B(nd,p0).  

   B2=zeros(n,1); 
   
   B2(1)        =random(t1,       1,1);   % the branching of the first generation 
   B2(2:B2(1)+1) =random(t1,      B2(1),1); % the branching of the second generation 
   B2(B2(1)+2:n)=binornd(nd,p0,n-B2(1)-1,1); % the branching of g=3 until g=ng with extinction 
 
 Binomial=B2; 
end 