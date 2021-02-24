%% Script for running minVpnorm
clear;
close all;
commandwindow;
clc;
% Domain size and discretization
domain = 'lbracket'; % 'biclamped'; % 'lbracket' % 
sizex = 200;
sizey = 200;
helem = 1.0;
% Optimization
penal = 3;
rmin = 5.5;
ft = 3;
Sytarget = 0.5;
comptarget = 1000;
betamax = 8;
betamaxsimul = 8;

% % Direct run
% solver.type = 1; 
% solver.cgtol1 = 1e-6;
% solver.cgtol2 = solver.cgtol1;
% filename = 'lbracket_direct_verify.mat';
% minVpnorm(domain,sizex,sizey,helem,penal,rmin,ft,Sytarget,comptarget,...
%     betamax,betamaxsimul,solver,filename);

%% Iterative runs
solver.type = 4; % This is for the MGCG solver
% Run tests
solver.cgtol1 = 1e-6;
while solver.cgtol1 < 2e-1
    solver.cgtol2 = solver.cgtol1;
    solver.cgmax1 = 10000;
    solver.cgmax2 = 10000;
    filename = sprintf('lbracket_%3.1e_%3.1e_reproduce.mat',solver.cgtol1,solver.cgtol2);
    % minVpnorm reproduces examples from section 5.1 and 5.2, l-bracket and bi-clamped
    % minVpnorm_long reproduces examples from section 5.3, u-bracket
    % minVpnorm_adpt reproduces examples from section 5.4, adaptive scheme
    minVpnorm(domain,sizex,sizey,helem,penal,rmin,ft,Sytarget,comptarget,...
        betamax,betamaxsimul,solver,filename);
    solver.cgtol1 = solver.cgtol1*10;
end
