function [ak_lag, data_prewhite, ak_ss]= nanprewhite_AR(data,varargin)

% Compute the first lag autocorrelation coefficient to prewhite the data as an AR(Kmax)function
% The prewhitened data are computed only if ak_lag is ss at the alpha_ak
% confidence limit. If it is not the case, original data are returned.

% Input:
%       data (array of floats)= data. Must be 1-D
% optional inputs: varargin
%       alpha_ak (integer)= confidence level for the first lag autocorrelation. Default=95%
%
% output:
%       ak_lag (float)= fist lag autocorrelation coefficient
%       data_prewhite (array of floats): data after removing of the first
%           lag autorcorrelation if ak_lag is ss, original data otherwise.
%       ak_ss (integer)= statistical significance of the first lag
%                        autocorrelation:alpha_ak is ss at the alpha_ak
%                        level, zero otherwise


%  Martine Collaud Coen, MeteoSwiss, 9.2020

% define the ak parameters
% data=data-nanmedian(data);

% Some sanity checks first
if isa(data,'float')==0 || min(size(data))>1
    error('the input "data" of Kendall_var has to be a 1-D array of floats');
end

% check optional input arguments
if ~varg_proof(varargin, {'alpha_ak'},true)
    return
end
% Set values from user input, or use defaults
alpha_ak = varg_val(varargin, 'alpha_ak', 95);
%sanity check
if alpha_ak>100
    error('the confidence limit has to be lower than 100%');
end

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
    %compute the autocorrelation
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