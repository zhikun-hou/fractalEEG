clear;
clc;
close all;


Hurst = 0.2;

[X,T] = fbm1d(Hurst,999,1);

% If fbm1d not work, try get it from "Fractional Brownian motion generator"
    % in matlab's FileExchange platform.

% =========================================================================

% test =====

% D_pred = tfd_Dubuc(struct(epsilon=10:100,visualize=true),X);

% D_pred = tfd_Higuchi(struct(epsilon=2.*[1:200]+1,visualize=true),X);

% D_pred = tfd_FA(struct(epsilon=1:500,visualize=true),X);

% D_pred = tfd_DFA(struct(epsilon=10:100,visualize=true),X);

D_pred = tfd_AFA(struct(epsilon=2.*[1:500]+1,visualize=true),X);


D_ground = 2-Hurst';
