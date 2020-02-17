function [ ] = Save_adjacency2file(file_number, nh,node,adj)
% This fucniton save the adjacency matrix. 
   adj1=zeros(node*node,1); 
   adj2=zeros(node*node+2,1);
   
   s1=node*node;
   adj1=reshape(adj,[s1,1]);
   adj2=[node;nh;adj1]; 
   filename = sprintf('%d%s',file_number,'.tree'); 
   dlmwrite(filename,adj2); 
end 