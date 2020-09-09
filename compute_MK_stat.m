function [result, S, vari,Z]=compute_MK_stat(time,data,resolution, varargin)
% compute all the components for the MK statistic
% in: time
% data
%resolution to compute the ties
% varargin:
        % alpha_MK= confidence level for the Mann-Kendall test in% 


%out: t= number of data per tie
% S= S statistic
% n=
%vari= kendall variance
%Z: Z statistic of the MK test
%P: probability for the statistical significance
%slope= Sen's slope
%slope_min=LCL
%slope_max=UCL

% check arguments
if ~varg_proof(varargin, {'alpha_MK','alpha_CL'},true)
    return
end

% Set values from user input, or use defaults
alpha_MK = varg_val(varargin, 'alpha_MK', 95);
alpha_CL = varg_val(varargin, 'alpha_CL', 90);

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

