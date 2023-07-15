% a script to run the script_CQ_timestepping to make convergence plots
% last modified: May 20, 2016

% don't forget to modify boundary conditions and lag time in the other
% script

external = 1;

k = 2; % polynomial degree
TT= 2; % final time

maxlev=7;
MM = 20*(1:maxlev);
levels = (1:maxlev);
errors = [];

for lev = levels
    M=MM(lev);
    script_CQ_timestepping_elasticity;
    disp(errors);
end

%%

loglog(levels,errors(:,1),'-o',...
       levels,errors(:,2),'-o',...
       levels,levels.^(-k),...
       levels,levels.^(-k-1))
xlabel('1/h')
ylabel('error')
legend('L2error','H1error','O(h^{k})','O(h^{(k+1)})',...
            'Location','northeast')
title('CQelasticty')