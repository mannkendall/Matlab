% Copyright 2020 MeteoSwiss, contributors of the original matlab version of the code listed in ORIGINAL_AUTHORS
% 
% Distributed under the terms of the BSD 3-Clause License.
% 
% SPDX-License-Identifier: BSD-3-Clause

function [slope,LCL,UCL]=Sen_slope(time,data,vari, varargin)

% Compute the Sen's slope, that is the median of the slopes for each
%   intervalles (xj-xi)/(j-i), j>i 
% The confidence limits are compute with an interpolation, which is very
%    important if the number of valid data is small (<=10).

%input: 
%       time (array of float)= time in datenum format
%       data (array of float)=  data. Must be 1-D array
%       vari (float) = Kendall variance (calculated with Kendall_var)
% Optional input: varargin: 
%       alpha_CL (float)= confidence level for the confidence limit in %, 
%           usually takes as90 or 95%. Default: 90% 

%output: 
%       slope (float)= Sen's slope
%       LCL (float)= lower confidence limit
%       UCL (float)= upper confidence limit

% Source:
%       Gilbert, 1987

% Some sanity checks first
if isa(time,'float')==0 || min(size(time))>1
    error('the input "time" of Sen_slope has to be a 1-D array of floats');
end
if isa(data,'float')==0 || min(size(data))>1
    error('the input "data" of Sen_slope has to be a 1-D array of floats');
end
if isa(vari,'float')==0 || max(size(vari))>1
    error('the input "vari" of Sen_slope has to be a float');
end

% check optional arguments
if ~varg_proof(varargin, {'alpha_CL'},true)
    return
end
% Set values from user input, or use defaults
alpha_CL = varg_val(varargin, 'alpha_CL', 90);
%sanity check:
if alpha_CL>100
    error('the confidence limit has to be lower than 100%');
end

%length of the time series
L=length(data);
% preallocation for the slopes between pairs
dij=NaN(L,L);

%compute the slope from all the pairs
for i=1:L-1
    dij(i+1:L,i)=(data(i+1:L)-data(i))./(time(i+1:L)-time(i));
end

%take the median
%put all the slopes in a single array
d=reshape(dij,1,L^2);
%remove NaN
d=d(~isnan(d));
N=size(d,2);
%compute the Sen's slope
slope=median(d);
% sort the slopes to compute the confidence limits
d_ascending=sort(d);

%compute the confidence limits for the defined confidence level alpha_CL
Cconf=-norminv((1-alpha_CL/100)/2)*vari^0.5;
M1=real((N-Cconf)/2);
M2=real((N+Cconf)/2);
if M1>1
    %LCL=d_ascending(round(real(M1)));
    %interpolation to obtain the best CL
    a1=d_ascending(ceil(M1));
    a2=d_ascending(floor(M1));
    LCL=a2+(a1-a2)*(M1-fix(M1));
elseif ~isempty(d_ascending)
    LCL=d_ascending(1); % the LCL cannot be smaller than the smallest slope
else
    LCL=NaN;
end
if M2>1
    %UCL=d_ascending(round(real(M2)));
    %interpolation to obtain the best CL
    a1=d_ascending(ceil(M2));
    a2=d_ascending(floor(M2));
    UCL=a2+(a1-a2)*(M2-fix(M2));
elseif ~isempty(d_ascending)
    UCL=d_ascending(end); % the LCL cannot be smaller than the smallest slop
else
    UCL=NaN; % the UCL cannot be larger than the largest slope
end

fclose('all');
