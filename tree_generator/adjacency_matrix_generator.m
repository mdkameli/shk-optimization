function [ adjt,nht,St,n1t,nodet] = adjacency_matrix_generator(B,ng)
node=0;nh5=0; S=zeros(1,ng); n1=zeros(1,ng);
 
   S(1) = B(1); n1(1)= S(1)+1;          
   S(2) = sum(B(2:n1(1))); n1(2)= S(2)+n1(1);     
   for i=1:ng-2
       S(i+2)=sum( B((n1(i)+1):(n1(i+1)))); 
       n1(i+2)=S(i+2)+n1(i+1);    
   end
   node=n1(ng); adj=zeros(node,node);f=0;
   for i=1:node
       if (1+B(i)+f) > node
        break
       end 
   for j=2+f:1+f+B(i)
       adj(i,j)=1;
   end  
        f=f+B(i);
   end 
 
   for i=1:node    
     for j=1:node
         adj(j,i)=adj(i,j);
     end 
   end

    %% find the heminodes (H)
   for i=2:node   % start from 2 to not confuse the central node as heminode in case it has only one branching
       if (sum(adj(i,1:node))==1)
         nh5=nh5+1; 
         bb(nh5)=i; % this matrix has the dimention of the peripherals and it does not matter where they are located. 
       end
   end  
adjt=adj; nht=nh5; St=S; n1t=n1; nodet=node; 
end 
