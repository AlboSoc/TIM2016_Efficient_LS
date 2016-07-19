%    LS_3PF   Evaluates 3-parameter LS fitting. Elements in matrix D0_T * D0
%                      are calculated directly based on the closed form
%                      of sums
%
%       S = LS_3PF(x, f_per_fs)
%
%       Input arguments:
%         x              - Input data
%                       
%         f_per_fs:     - relative signal frequency
%                       (equals to f_signal / f_sampling) 
        
%       
%       Output arguments:
%         S -    Structure containing the parameters that are estimated,
%                that is, cosine and sine amplitudes and offset (A, B and
%                C), respectively

% Written by Balázs Renczes
% $Id: phase_eval_prec.m,v 1.0 20-May-2016 $
% Copyright (c) 2016 by Balázs Renczes
% All rights reserved.

function S = LS_3PF(x, f_per_fs)

N = length(x);
phi1 = 2*pi*f_per_fs;

if size(x,2) > 1,  x = x';  end

% H = D0_T * D0;

% Constructing the elements of H

h11 = N/2 + cos((N-1)*phi1) * sin(N*phi1) / (2*sin(phi1));

h12 = 0.5 * sin(N*phi1) * sin((N-1)*phi1) / sin(phi1);

h13 = cos((N-1)*phi1/2)* sin(N*phi1/2)/ sin(phi1/2);

h22 = N - h11;

h23 = sin((N-1)*phi1/2)* sin(N*phi1/2)/ sin(phi1/2);

h33 = N;

H = [h11 h12 h13; h12 h22 h23; h13 h23 h33];

%Constructing the elements of D0

D0 =  [cos(phi1*(0:(N-1)))'  sin(phi1*(0:(N-1)))'  ones(N,1)];

s = H^-1 * (D0' * x);

S.A = s(1);
S.B = s(2);
S.C = s(3);
