function [Ud] = postComputeDirichlet3D(tspan,T,k,g0)
%
% [Ud] = postComputeDirichlet3D(tspan,T,k,g0)         
%
% Input:
% tspan          : vector of timesteps
% T              : Data structure. Basic FEM tetrahedrization, updated
% k              : scalar. Polynomial degree
% g0             : array of vectorized functions  ( of (t,x,y,z) )
%
% Output:
% Ud            : N_g0(N_dir) x #timespan matrix of Dirichlet test values
% Last modified July 8, 2015.

[~,~,~,~,~,dir,~]=bdDOF3D(T,k);
nT=length(tspan);
nD=length(dir);
for i=1:3
    for tk=1:nT
        t=tspan(tk);
        g{i}{tk}= @(x,y,z) g0{i}(t,x,y,z);
    end
end
parfor i=1:3
   Udcomp{i}=dirichletBC3D(g{i},T,k)
   Udcomp{i}=Udcomp{i}(dir,:);
end
Ud=[Udcomp{1};Udcomp{2};Udcomp{3}];