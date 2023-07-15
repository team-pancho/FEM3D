% a script to run the script_CQ_timestepping to make convergence plots
% last modified: May 9, 2016

% don't forget to modify boundary conditions and lag time in the other
% script

external = 1;

k = 2;
TT= 2;

MM = 5*2.^(1:8);
levels = (1:8);
errors = [0 0];

for lev = levels
    M=MM(lev);
    script_CQ_timestepping;
    disp(errors);
end