function dataPW=prewhite_D(data, colonne, resolution,varargin)

%calculate all the necessary prewhitened datasets to asses the statistical
%significance and to compute the slope the slope

%in: data as timetable, matrix or structure
% colonne: to be suppressed
%resolution= is used to compute the number of ties
%varargin: alpha_ak = statistical significance on % for the first lag
%                      autocorrelation. Default=95

%out: timetable with 3 PW timeseries:
%       data.PW = PW with the first lag autocorrelation of the data
%       data.PW_cor=PW corrected with 1/(1-ak1)
%       data.TFPW_WS = PW with the first lag autocorrelation of the data after
%                   detrending computed from PW data (see Wang & Swail)
%       data.TFPW_Y= method of Yue et al 2002, not 1/1-ak1) correction,
%                       detrend on original data
%       data.VCTFPW = PW with the first lag autocorrelation of the data
%                   after detrending + correction of the PW data for the variance (see
%                   Wang 2015)

%Martine Collaud Coen, 9.2020

% check arguments
if ~varg_proof(varargin, {'alpha_ak'},true)
    return
end

% Set values from user input, or use defaults
alpha_ak = varg_val(varargin, 'alpha_ak', 95);

%determination of the data to analyse
if isstruct(data)
    time=data.start_time;
    data=data.colonne;
elseif istimetable(data)
    time=datenum(data.Time);
    data=data.(colonne);
else
    time=datenum(data(:,1:6));
    data=data(:,colonne);
end
data(abs(data)==inf)=NaN;
nb_loop=0;
%create the timetable
dataPW=timetable(datetime(datevec(time)));
%calcul de l'autocorrelation:
[c.PW, dataARremoved, c.ss]= nanprewhite_ARok(data,'alpha_ak',alpha_ak);


% compute data PW corrected
if sum(~isnan(dataARremoved))>0 & c.ss==alpha_ak & c.PW>=0.05
    dataPW.PW= dataARremoved;
    dataPW.PW_cor= dataARremoved./(1-c.PW);
    
    % data VCTFPW corrected
    %calcul of the trend slope o the PW data
    
    t=Nb_tie_D(dataARremoved./(1-c.PW),resolution);
    [~,n]=S_test(dataARremoved./(1-c.PW), time);
    vari=Kendall_var(dataARremoved./(1-c.PW),t,n);
    [b0PW,~,~]=Sen_slope(time,dataARremoved./(1-c.PW),vari); %slope of PW data
    
    t=Nb_tie_D(data,resolution);
    [~,n]=S_test(data, time);
    vari=Kendall_var(data,t,n);
    [b0or,~,~]=Sen_slope(time,data,vari); %slope of original data
    
    %remove the trend
    dataDetrendPW=data-b0PW*(time-time(1));
    dataDetrendor=data-b0or*(time-time(1));
    %compute the autocorrelation of the de-trended time series:
    [c.VCTFPW, dataARremovedor, c.ssVC]= nanprewhite_ARok(dataDetrendor,'alpha_ak',alpha_ak);
    
    [akPW, dataARremovedPW, ssPW]= nanprewhite_ARok(dataDetrendPW,'alpha_ak',alpha_ak);
    
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
        
        t=Nb_tie_D(dataARremovedPW,resolution);
        [~,n]=S_test(dataARremovedPW, time);
        vari=Kendall_var(dataARremovedPW,t,n);
        [b1PW,~,~]=Sen_slope(time,dataARremovedPW,vari);
        
        %remove the trend
        
        while  abs(akPW-c1)>0.0001 & abs(b1PW-b0PW)>0.0001
            if akPW>=0.05 & ssPW==95
                nb_loop=nb_loop+1;
                
                dataDetrendPW=data-b1PW*(time-time(1));
                c1=akPW; b0PW=b1PW;
                [akPW, dataARremoved2PW,ssPW]= nanprewhite_ARok(dataDetrendPW,'alpha_ak',alpha_ak);
                if akPW>0 & ssPW==95
                    dataARremoved2PW=[data(1);(data(2:end)-akPW*data(1:end-1))./(1-akPW)];
                    t=Nb_tie_D(dataARremoved2PW,resolution);
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
    %c.PW=0;
    dataPW.TFPW_Y= data;
     dataPW.TFPW_WS= data;
    %c.TFPW=0;
    dataPW.VCTFPW= data;
    %c.VCTFPW=0;
% % %     c.VCTFPW=NaN;
% % %     c.ssVC=NaN;
% % %     c.TFPW_WS=NaN;
end

fclose('all');