function [ Tri ] = General_Binary_Branching(n,m0,m1,m2,p0)
 %General_Binary_Branching: The zero branching for the root is either 1 or 2, and for the rest, it is 0, 1, or 2. 
 % p(m1)=p0, p(m1)=p(m2)=(1-p0)/2
 B2=zeros(n,1); B1=rand(n,1); 
  if (B1(1) <= 0.5)
      B2(1)=m1; 
  else
      B2(1)=m2; 
  end 
 
 for i=2:n
     if B1(i) <p0 
         B2(i)=m0; 
     elseif (B1(i) >= p0) && (B1(i) < (p0+1)/2)
         B2(i)=m1; 
     else
         B2(i)=m2; 
     end 
 end
 Tri=B2; 
end 