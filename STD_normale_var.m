function Z=STD_normale_var(S,var_S)
% calculate the normal standard variable  Z
%see Gilbert 1987

%in: S is the S statistic of the Mann-Kendall test computed from S_test
%      it is a entire reprensenting the difference between the position and
%      negative pais in the times series (TO DO: reformuler)
% var-S if the variance of the time series taking into account the ties in
%       values and time. It is computed by Kendall_var

%out: Z is the S statistic weighted by the variance
%Martine Collaud Coen, february 2019
if S==0
    Z=0;
elseif S > 0
    Z=(S-1)/((var_S)^0.5);
else
    Z=(S+1)/((var_S)^0.5);
end
 fclose('all');