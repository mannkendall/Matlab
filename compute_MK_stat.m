function [result, S, vari,Z]=compute_MK_stat(time,data,resolution, varargin)
% compute all the components for the MK statistic of the main MK function

% input:
%       time (array)= time in datetime format
%       data (array of floats)= data. Must be 1-D
%       resolution (float)= interval needed to compute the ties (Nb_tie_D
% Optional input: varargin:
%       alpha_MK (float)= confidence level for the Mann-Kendall test in %. Default = 95%
%       alpha_CL (float)= confidence level for the confidence levels of the Sen's slope
%       in %. Default = 90%

%output:
%       result is a structure comprising:
%           result.P (float)= probability for the MK test
%           result.ss (float)= statistical significance of the MK test.
%               Is equal to alpha_MK if ss and zero if not ss.
%           result.slope (float)= Sen's slope
%           result.UCL (float)= upper confidence limit of the Sen's slope
%           result.LCL (float)= lower confidence limit of the Sen's slope
%       S (integer)= S statistic computed from S_test
%       vari (array of integer)= Kendall variance computed from Kendall_var
%       Z (float)= S statistic weighted by the Kendall variance, computed
%           from STD_normal_var

% Some sanity checks first
if isdatetime(time)==0 || min(size(time))>1
    error('the input "time" of compute_MK_stat has to be a 1-D array of datetime ');
end
if isa(data,'float')==0 || min(size(data))>1
    error('the input "data" of compute_MK_stat has to be a 1-D array of floats');
end
if isa(resolution,'float')==0 || max(size(resolution))>1
    error('the input "resolution" of compute_MK_stat has to be a single float');
end

% check optional arguments
if ~varg_proof(varargin, {'alpha_MK','alpha_CL'},true)
    return
end
% Set values from user input, or use defaults
alpha_MK = varg_val(varargin, 'alpha_MK', 95);
alpha_CL = varg_val(varargin, 'alpha_CL', 90);

% sanity check
if alpha_MK>100 || alpha_CL>100
    error('confidence limit should be smaller then 100%');
end

t=Nb_tie_D(data,resolution);
[S,n]=S_test(data,datenum(time));
vari=Kendall_var(data,t,n);
Z=STD_normale_var(S,vari) ;
if sum(~isnan(data)) > 10
    result.P=2*(1-normcdf(abs(Z),0,1));
else
    load ('Prob_MK_n');
    result.P=Prob_MK_n(abs(S)+1,sum(~isnan(data)));
end
%determine the ss
if result.P<=1-alpha_MK/100
    result.ss= alpha_MK;
else
    result.ss= 0;
end
[slope,slope_min,slope_max]=Sen_slope(datenum(time),data,vari,'alpha_CL',alpha_CL);
result.slope=slope.*365.25;
result.UCL=(slope_max.*365.25)';
result.LCL=(slope_min.*365.25)';

