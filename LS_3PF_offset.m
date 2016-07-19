%    LS_3PF_offset   Evaluates 3-parameter LS fitting with time offsetting.
%                    Elements in matrix D0_T * D0
%                      are calculated directly based on the closed form
%                      of sums
%
%       S = LS_3PF_offset(x, f_per_fs)
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
% $Id: LS_3PF_offset.m,v 1.0 20-May-2016 $
% Copyright (c) 2016 by Balázs Renczes
% All rights reserved.

function S = LS_3PF_offset(x, f_per_fs)

N = length(x);
phi1 = 2*pi*f_per_fs;

if size(x,2) > 1,  x = x';  end


% H = D0_T * D0;

% Constructing the elements of H

h11 = N/2 + 0.5* sin(N*phi1) / sin(phi1);

h12 = 0;

h13 = sin(N*phi1/2)/ sin(phi1/2);

h22 = N - h11;

h23 = 0;

h33 = N;

H = [h11 h12 h13; h12 h22 h23; h13 h23 h33];

%In order to ensure optimal conditioning, the last row and column of H
% has to be divided by sqrt(2)

H(:,3) = H(:,3) / sqrt(2);
H(3,:) = H(3,:) / sqrt(2);

%Constructing the elements of D0

l = (N-1)/2;   %time offsetting

D0 =  [cos(phi1*((0:(N-1))-l))'  sin(phi1*((0:(N-1))-l))'  ones(N,1)];

D0(:,3) = D0(:,3) / sqrt(2);

s = H^-1 * (D0' * x);

%Calculating the original parameters

S.A = s(1)*cos(2*pi*f_per_fs*l) - s(2)*sin(2*pi*f_per_fs*l);
S.B = s(1)*sin(2*pi*f_per_fs*l) + s(2)*cos(2*pi*f_per_fs*l);
S.C = s(3)/sqrt(2);
