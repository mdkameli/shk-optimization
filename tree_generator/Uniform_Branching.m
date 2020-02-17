function [ Uniform ] = Uniform_Branching(n,nd)
%Uniform_Branching: Zero branching is allowed after the second generation thus the branching at the root is (1-nd). 
%For the rest, it is (0-nd). The branching is generated using uniform random numbers. 
   B2=zeros(n,1); 

   B2(1)        =randi([1 nd],        1,1);   % the branching of the first generation 
   B2(2:B2(1)+1) =randi([1 nd],     B2(1),1); % the branching of the second generation 
   B2(B2(1)+2:n)=randi([0 nd],n-B2(1)-1,1); % the branching of g=3 until g=ng 
 
 Uniform=B2; 
end 