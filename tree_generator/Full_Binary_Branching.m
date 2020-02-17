function [ Bin ] = Full_Binary_Branching(n,m0,m2,p0)
%Full_Binary_Branching: The branching is m2=2 until generation 2, and it is either 0 or 2.
% The probability of zero branching is p0. 

    B=zeros(n,1);     
    B(1:m2+1)=m2; 
     for k5=B(1)+2:n
        z1=rand(); 
        if z1 <p0 
            B(k5)=m0; 
        else
            B(k5)=m2;
        end  
    end 
 Bin=B; 
end 