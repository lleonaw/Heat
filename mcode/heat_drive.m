%%  Driver for heat.m 
%%  Heat equation, no preferred direction, central flux
%  This code is for discontinuous Galerkin method
%  Explicit scheme, SSP-RK2

%   clear all; format long;

close; 
global Ne Nx
Ne = 20;          % Number of elem
N = 4;            % Poly. order
Nx = N + 1;       % Numb of points in each elem.

succ = heat; 

