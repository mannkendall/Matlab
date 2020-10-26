
% Copyright 2020 MeteoSwiss, contributors of the original matlab version of the code listed in ORIGINAL_AUTHORS
% 
% Distributed under the terms of the BSD 3-Clause License.
% 
% SPDX-License-Identifier: BSD-3-Clause

function dataPW=prewhite(data, colonne, resolution,varargin)

%calculate all the necessary prewhitened datasets to asses the statistical
%significance and to compute the Sen's slope for each of the prewhitening
%method, including 3PW

% inputs: 
%       data (timetable with only one argument): data to analyse
%       ?? colonne (integer or string) : define the data to use, TODO:to be suppressed
%       resolution (float)= the measurement resolution, i.e. delta value below which 2 measurements
%            are considered equivalent. It is used to compute the number of ties
%
% optional inputs: varargin: 
%       alpha_ak = statistical significance in % for the first lag
%                  autocorrelation. Default=95

% output: 
%       dataPW (timetable): contains 5 PW timeseries:
%           data.PW = PW with the first lag autocorrelation of the data
%           data.PW_cor=PW corrected with 1/(1-ak1)
%           data.TFPW_WS = PW with the first lag autocorrelation of the data after
%                   detrending computed from PW data (see Wang & Swail)
%           data.TFPW_Y= method of Yue et al 2002, not 1/1-ak1) correction,
%                       detrend on original data
%           data.VCTFPW = PW with the first lag autocorrelation of the data
%                   after detrending + correction of the PW data for the variance (see
%                   Wang 2015)

%Martine Collaud Coen, 9.2020

% First some sanity checks
if isa(data,'timetable')==0 || min(size(data))>1
    error('the input "data" of prewhite_D has to be a timetable');
end
if isa(resolution,'float')==0 || max(size(resolution))>1
    error('the input "resolution" of Nb_tie_D has to be a single float');
end
 
% check optional arguments
if ~varg_proof(varargin, {'alpha_ak'},true)
    return
end
% Set values from user input, or use defaults
alpha_ak = varg_val(varargin, 'alpha_ak', 95);
if alpha_ak>100
    error('the confidence limit has to be lower than 100%');
end

%determination of the array of time and of data
 time=datenum(data.Time);
    data=data.(colonne);
%check for the presence of infinites
data(abs(data)==inf)=NaN;
% counter of the number of loops
nb_loop=0;
%create the output timetable
dataPW=timetable(datetime(datevec(time)));
%calcul de l'autocorrelation:
[c.PW, dataARremoved, c.ss]= nanprewhite_AR(data,'alpha_ak',alpha_ak);

% compute data PW corrected
if sum(~isnan(dataARremoved))>0 & c.ss==alpha_ak & c.PW>=0.05
    dataPW.PW= dataARremoved;
    dataPW.PW_cor= dataARremoved./(1-c.PW);
    
    % data VCTFPW corrected
    %calcul of the trend slope of the PW data
    
    t=Nb_tie(dataARremoved./(1-c.PW),resolution);
    [~,n]=S_test(dataARremoved./(1-c.PW), time);
    vari=Kendall_var(dataARremoved./(1-c.PW),t,n);
    [b0PW,~,~]=Sen_slope(time,dataARremoved./(1-c.PW),vari); %slope of PW data
    
    t=Nb_tie(data,resolution);
    [~,n]=S_test(data, time);
    vari=Kendall_var(data,t,n);
    [b0or,~,~]=Sen_slope(time,data,vari); %slope of original data
    
    %remove the trend
    dataDetrendPW=data-b0PW*(time-time(1));
    dataDetrendor=data-b0or*(time-time(1));
    %compute the autocorrelation of the de-trended time series:
    [c.VCTFPW, dataARremovedor, c.ssVC]= nanprewhite_AR(dataDetrendor,'alpha_ak',alpha_ak);
    
    [akPW, dataARremovedPW, ssPW]= nanprewhite_AR(dataDetrendPW,'alpha_ak',alpha_ak);
    
    %calcul of TFPW correction Yue et al., 2002
      %blended data
    if sum(~isnan(dataARremovedor))>0
        dataPW.TFPW_Y= dataARremovedor+b0or*(time-time(1));
    else
        dataPW.TFPW_Y= data;
    end
    
    %calcul of TFPW correction Wang and Swail
    
    %avant C1>=0.05
    
    if abs(akPW)>=0.05 & ssPW==95
        %change so that while can be used with the same variable: ak is new c and c1 is last c
        %ak=c.VCTFPW; 
        c1=c.PW;
        
        dataARremovedPW=[data(1); ((data(2:end)-akPW*data(1:end-1))./(1-akPW))];
        
        t=Nb_tie(dataARremovedPW,resolution);
        [~,n]=S_test(dataARremovedPW, time);
        vari=Kendall_var(dataARremovedPW,t,n);
        [b1PW,~,~]=Sen_slope(time,dataARremovedPW,vari);
        
        %remove the trend
        
        while  abs(akPW-c1)>0.0001 & abs(b1PW-b0PW)>0.0001
            if akPW>=0.05 & ssPW==95
                nb_loop=nb_loop+1;
                
                dataDetrendPW=data-b1PW*(time-time(1));
                c1=akPW; b0PW=b1PW;
                [akPW, dataARremoved2PW,ssPW]= nanprewhite_AR(dataDetrendPW,'alpha_ak',alpha_ak);
                if akPW>0 & ssPW==95
                    dataARremoved2PW=[data(1);(data(2:end)-akPW*data(1:end-1))./(1-akPW)];
                    t=Nb_tie(dataARremoved2PW,resolution);
                    [~,n]=S_test(dataARremoved2PW, time);
                    vari=Kendall_var(dataARremoved2PW,t,n);
                    [b1PW,~,~]=Sen_slope(time,dataARremoved2PW,vari);
                    dataARremovedPW=dataARremoved2PW;
                    if nb_loop>10
                        break
                    end
                else
                    break
                end
            end
        end
    else
        b1PW=b0PW; akPW=c.PW;
    end
    %blended data
    
    if sum(~isnan(dataARremovedPW))>0
        dataPW.TFPW_WS= dataARremovedPW;
        c.TFPW_WS=akPW;
    else
        dataPW.TFPW_WS= data;
    end
    
    %correction VCTFPW
     % correction of the variance
    var_data=nanvar(data);
    var_data_TFPW=nanvar(dataARremovedor);
    dataARremoved_var=dataARremovedor.*(var_data/var_data_TFPW);
    %modify slope estimator (correction of the slope for the
    %autocorrelation
    if c.VCTFPW>=0
        bVC=b0or/sqrt((1+c.VCTFPW)/(1-c.VCTFPW));
        
    elseif c.VCTFPW<0
        bVC=b0or;
    end
    %add the trend again
    dataPW.VCTFPW= dataARremoved_var+bVC*(time-time(1));
    
else %no s.s. autocorrelation
    dataPW.PW= data;
    dataPW.PW_cor= data;
    dataPW.TFPW_Y= data;
     dataPW.TFPW_WS= data;
    dataPW.VCTFPW= data;
end

fclose('all');