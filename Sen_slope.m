function [slope,LCL,UCL]=Sen_slope(time,data,vari, varargin)

%Sen's slope, that is the median of the slopes for each
%intervalles (xj-xi)/(j-i), j>i (see Gilbert 1987)
% The confidence limits are compute with an interpolation
% that is important is the number of data is small, such as yearly averages for a 10 year�s trend.

%in: time= time in datenum
%   data= one column of data
%   vari = Kendall variance (calculated with Kendall_var)
%   varargin: alpha_CL= confidence level for the confidence limit (90 or 95%, if
%   not specified, 90%)

%out: slope= median==Sen slope
%     LCL= lower confidence limit
%     UCL= upper confidence limit

% Martine Collaud Coen, february 2019

% check arguments
if ~varg_proof(varargin, {'alpha_CL'},true)
    return
end

% Set values from user input, or use defaults
alpha_CL = varg_val(varargin, 'alpha_CL', 90);

L=length(data);

dij=NaN(L,L);

%compute the slope from all the pairs
for i=1:L-1
    dij(i+1:L,i)=(data(i+1:L)-data(i))./(time(i+1:L)-time(i));
end

%take the median
d=reshape(dij,1,L^2);
d=d(~isnan(d));
N=size(d,2);
slope=median(d);
d_ascending=sort(d);

%compute the confidence  limits for the defined confidence level alpha_CL
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
