function Z=STD_normale_var(S,var_S)
% calculate the normal standard variable  Z

%input: 
%       S (integer)= S statistic of the Mann-Kendall test computed from
%           S_test 
%       var_S (float)= variance of the time series taking into account the ties in
%           values and time. It is computed by Kendall_var

%output: 
%       Z (float)= S statistic weighted by the variance

% Source:
%       Gilbert 1987

% sanity checks first
if isa(S,'float')==0 || max(size(S))>1
    error('the input "S" of STD_normale_var has to be a 1-D float');
end
if isa(var_S,'float')==0 || max(size(var_S))>1
    error('the input "var_S" of STD_normale_var has to be a 1-D float');
end

%compute Z
if S==0
    Z=0;
elseif S > 0
    Z=(S-1)/((var_S)^0.5);
else
    Z=(S+1)/((var_S)^0.5);
end
 fclose('all');