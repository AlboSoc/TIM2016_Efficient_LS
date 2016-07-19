%    LS_4PF_offset   Evaluates 4-parameter LS fitting with time offsetting.
%                    Elements in matrix D0_T * D0
%                      are calculated directly based on the closed form
%                      of sums
%
%       S = LS_4PF_offset(x, f_per_fs_est)
%
%       Input arguments:
%         x              - Input data
%                       
%         f_per_fs_est:     - relative signal frequency estimator  (optional)
%                       (equals to f_signal / f_sampling) 
        
%       
%       Output arguments:
%         S -    Structure containing the parameters that are estimated,
%                that is, cosine and sine amplitudes,offset and relative
%                signal frequency (A, B, C and f_per_fs), respectively

% Written by Balázs Renczes
% $Id: LS_4PF_offset.m,v 1.0 20-May-2016 $
% Copyright (c) 2016 by Balázs Renczes
% All rights reserved.

function S = LS_4PF_offset(x, f_per_fs)

N = length(x);
if size(x,2) > 1,  x = x';  end

if nargin <2
    
    X = abs(fft(x));
    [dummy, f_per_fs] = max(X(2:end));   %calculation of frequency by FFT
    f_per_fs = f_per_fs/N;
end

% Initial parameter estimation
phi1 = 2*pi*f_per_fs;

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

l = N/2 -1;   %time offsetting

D0 =  [cos(phi1*((0:(N-1))-l))'  sin(phi1*((0:(N-1))-l))'  ones(N,1)/sqrt(2)];

s = H^-1 * (D0' * x);

A_comma = s(1);
B_comma = s(2);
C_comma = s(3) / sqrt(2);

n_iter = 0;

while 1
 
    n_iter = n_iter + 1;
    
  % Constructing the elements of H  
phi1 = 2*pi*f_per_fs;
  
h11 = N/2 + 0.5* sin(N*phi1) / sin(phi1);

h12 = 0;

h13 = sin(N*phi1/2)/ sin(phi1/2);

h14 = -A_comma/4 * (-N*cos(N*phi1)/sin(phi1) + sin(N*phi1) * cos(phi1)/(sin(phi1))^2 );

h22 = N - h11;

h23 = 0;

h24 =  B_comma/4 * (-N*cos(N*phi1)/sin(phi1) + sin(N*phi1) * cos(phi1)/(sin(phi1))^2 );

h33 = N;

h34 = -A_comma/2 * (-N * cos(N*phi1/2)/sin(phi1/2) + sin(N*phi1/2)*cos(phi1/2)/(sin(phi1/2))^2);

h44 = (A_comma^2+B_comma^2)*(N^3-N)/24 + (B_comma^2-A_comma^2) * (N^2/8 * sin(N*phi1)/sin(phi1) + 1/4 * (sin(N*phi1)/2* (-1/sin(phi1)-2*cos(phi1)^2/sin(phi1)^3))+...
    N/4* cos(N*phi1) * cos(phi1)/ sin(phi1)^2  );

H = [h11 h12 h13 h14; h12 h22 h23 h24; h13 h23 h33 h34; h14 h24 h34 h44];

gamma = sqrt((A_comma^2+B_comma^2)*(N^3-N)/24);

H(:,3) = H(:,3) / sqrt(2);
H(3,:) = H(3,:) / sqrt(2);

H(:,4) = H(:,4)/gamma;
H(4,:) = H(4,:)/gamma;
  
D =  [cos(phi1*((0:(N-1))-l))'  sin(phi1*((0:(N-1))-l))'  ones(N,1)/sqrt(2)  (-A_comma*((0:(N-1))-l)'.*sin(phi1*((0:(N-1))-l))'+ B_comma*((0:(N-1))-l)'.*cos(phi1*((0:(N-1))-l))')/gamma];

s = H^-1 * (D' * x);

A_comma = s(1);
B_comma = s(2);
C_comma = s(3)/sqrt(2);

f_per_fs = f_per_fs + s(4)/gamma/(2*pi);

if abs(s(4)/gamma/(2*pi) / f_per_fs) < 1e-7  break; end
if n_iter >= 6 break;  end

end

err = x - (A_comma*cos(2*pi*f_per_fs*((0:(N-1))-l)) + B_comma*sin(2*pi*f_per_fs*((0:(N-1))-l)) + C_comma)';


%Calculating the original parameters

S.A = s(1)*cos(2*pi*f_per_fs*l) - s(2)*sin(2*pi*f_per_fs*l);
S.B = s(1)*sin(2*pi*f_per_fs*l) + s(2)*cos(2*pi*f_per_fs*l);
S.C = C_comma;
S.f = f_per_fs;
S.erms = 1/N* err'*err;
