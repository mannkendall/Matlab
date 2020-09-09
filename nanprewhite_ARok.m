function [ak_lag, data_prewhite, ak_ss]= nanprewhite_ARok(data,varargin)

%define the first lag autocorrelation coefficient to prewhite the data as an AR(Kmax)function
%
% Input:
% data: one column of data
%varargin
% alpha_ak= confidence level for the first lag autocorrelation in %
% fig: = 1 if the autocorrelogram has to be plotted, 0 for no figure
%
% output:
%
%ak_lag, partial coefficients of AR(Kmax) process
%K: degree of autocorrelation process
% data_prewhite: data after removing of the first lag autorcorrelation
%ak_ss= ss of the first lag autocorrelation


%  Martine Collaud Coen, MeteoSwiss, 9.2019

% define the ak parameters
%data=data-nanmedian(data);

% check arguments
if ~varg_proof(varargin, {'alpha_ak'},true)
    return
end

% Set values from user input, or use defaults
alpha_ak = varg_val(varargin, 'alpha_ak', 95);

%%% TO DO check if data is a one dimentional array of double.

if sum(isnan(data))==length(data)
    ak_lag=NaN; data_prewhite=NaN;
else
    %remove inf if existing
    data(abs(data)==Inf)=NaN;
    %specifies the number of lags supposed to have a significant AC coefficient
    p=5;
    % number of lag to be computed
    nblag=10;
    %restrict the number of computed lag and the number of significant AC coef. if the time series is short
    if nblag>length(data)/2
        nblag= floor( length(data)/2);
    end
    if p>=nblag
        p=nblag-1;
    end
    [x,~]=nanautocorr(data,nblag,p);
    
    % compute the  confidence limits for the autocorrelation
    x1 = double(x)./sum(~isnan(data));
    [~,~,K] = levinson(x1,p);
    ak.coef = -K;
    uconf=norminv(1-(1-alpha_ak/100)/2)/sqrt(sum(~isnan(data))); %lconf=-uconf;
    
     
    ak_lag=x(2);
    %compute the prewhitened data only is ak(1) is s.s. at 95% confidence limit 
    %other wise the orgininal is given as an output
    if abs(ak.coef(1))-uconf<0
        data_prewhite=data;
        ak_ss=0; 
        warning('no s.s. autocorrelation');
    else
        y=NaN(length(data),1);
        y(2:end,1)=x(2)*data(1:end-1);
        data_prewhite=(data-y(:,1));
        ak_ss=alpha_ak;
    end
end