% This program generate an ensnmble of random trees, the outputs are number
% of heminode, total nodes, adjacency matrices and ....
clear all; clc;
it=1000;                 % Number of trails
p0=0.5 ;                 % Probability of zero branching 
m0=0; m1=1; m2=2;        % The bracings for full and genreal binary trees
ng=4;                    % The maximum number of allowed generations
nd=2;                    % Mmaximum nmber of daughters 
nl= nd^ng;               % Maximum number of perhiphrals (Cayley tree)
nn=(nd^(ng+1)-1)/(nd-1); % Maximum total number of the nodes (Cayley tree)
 binomail1=makedist('Binomial',nd,p0); % make the distribution for Binomial p0, and N=nd
truncate_binomial=truncate(binomail1,1,inf);           % truncate(bi1,lower,upper), truncated binomial distrubtion (exlcude the zero branching). 

for ii=1:it
   %% Choose your favorite branching from the following:
 %B=Full_Binary_Branching(nn,m0,m2,p0);               % Full Binary trees, Branching is 0(m0) or 2(m2), (set ng=4)
 %B=General_Binary_Branching(nn,m0,m1,m2,p0);        % General Binary trees, Branching is 0 (m0),1(m1), or 2(m2), (set ng=3)
 %B=Uniform_Branching(nn,nd);                        % Uniform Branching, zero branching start after g>=2, nd=4,ng=4
 B=Binomial_Branching(nn,nd,truncate_binomial,p0);  % Binomail branching, zero brnaching is start after g>=2, nd=4,ng=4
%% Adjacency Matrix 
 [adj,nh,S,n1,node]=adjacency_matrix_generator(B,ng);     % constructing the adjacency matrix from the branching
 %Save_adjacency2file(ii, nh,node,adj);                  % save the trees as file_number.tree
%% Laplacian 
 deg=zeros(node,node);
 for j=1:node 
     deg(j,j)=sum(adj(j,:)); %Degrees of each node
 end 
 laplacian=deg-adj;          %Laplacian matrix
 eig1=eig(laplacian);
 lap1(ii)=eig1(2);           %Save second eigenvalue of the laplaxcian 
%% Saving tree properties to arryas
   r1=node-S(ng); 
   k2=0; aa=find(B(1:r1)==0); k2=length(aa);    
   
   nh_inside(ii)=k2;          %The number of heminodes that are in g<ng 
   nh_g4(ii)=S(ng);           %The number of node in g=ng, 
   T(ii,1:ng)=S(1:ng);        %The population of shells 1-ng
   T(ii,ng+1)=nh;             %The total number of heminodes
   T(ii,ng+2)=node;           %The total number of nodes
   T(ii,ng+3)=nh/node;        %The ratio R=nh/node
   T(ii,ng+4)=nh/(node*node); %The ratio nh/(node^2)
   N1(ii)=node; H1(ii)=nh; 
end
%%
clf
hist(N1,100);  
title('Histogram of the total number of the nodes')