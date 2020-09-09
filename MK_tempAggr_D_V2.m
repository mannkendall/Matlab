function result=MK_tempAggr_D_V2(data_tempAgg, PW_method, resolution, varargin)

% The MK test and the Sen slope is applied on the given time granularity
% Three prewhitening methods are applied. PW (Yue et al., 2002) and TFPW_Y(trend free PW, Wang and Swail, 2001)
%to compute the statistical significance and VCTFPW (Wang, W., Chen, Y., Becker, S., & Liu, B. (2015). Variance Correction Prewhitening Method for Trend Detection in Autocorrelated Data. J. Hydrol. Eng., 20(12), 4015033-1-04015033�10. https://doi.org/10.1061/(ASCE)HE.1943-5584.0001234.)
%to compute the Sen's slope
% only the statistically significant autocorrelation are taken into account
% for the prewhitening.
% The ss of the trends is taken at 95% confidence limit
% The upper and lower confidence limits are given by the 90% of the all intervals differences distribution.
% The significance level is given by the MK test and has therefore no direct
%relation to the confidences limits.
% if seasonal Mann-Kendall is applied, the yearly trend is assigned only if
% the results of the seasonal test are homogeneous. In case of not ss

%references:  WMO-GAW publication N. 133, annexe E, p. 26
% and the explanations of  MULTMK/PARTMK de C. Libiseller
% and the book ofof Gilbert 1998


%IN
% data_tempAgg is a timetable, or a structure of timetables if the seasonal
% Mann-Kendall test with a temporal segmentation is used. The timetable should have
% only one field (apart the time).
% PW_method is the used PW method (3PW; PW, TFPW_Y, TFPW_WS, VCTFPW). Default is 3PW.
% resolution is taken into account to determine the number of ties. I
% always take it similar to the resolution of the instrument or a little
% bit higher. This parameters can be determinant for the results but not very
% sensitive.
% varargin: alpha_MK: confidence limit for Mk test in %. Default value is 95%
%           alpha_CL: confidence limit for the confidence limits of the Sen's slope in %. Default value is 90%
%           alpha_Xhomo: confidence limit for the homogeneity between seasons in %. Default value is 90%
%           alpha_ak: confidence limit for the first lag autocorrelation in %. Default value is 95%

%OUT
%The result is a table  with the following fields
%result.period:year or 12 monthes + year or 4 meteorological seasons+year (TO DO ????)
%result.P: probability for the statistical significance. If 3PW is applied,
%P= max(P_PW, P_TFPW_Y);
%result.ss= statistical significance: alpha% if the test is ss at the alpha
%confidence level. Default=95%
%                                     0 if the test is not ss at the alpha
%                                       confidence level
%                                     -1 if the test is a TFPW_Y false
%                                       positive at alpha% confidence level
%                                     -2 if the test is a PW false
%                                       positive at alpha% confidence level
%result.slope: Sen's slope in units/y
%result.UCL: upper confidence level in units/y
%result.LCL: lower confidence level in units/

% Martine Collaud Coen , MeteoSwiss, 9.2020

% check arguments
if ~varg_proof(varargin, {'alpha_MK','alpha_CL','alpha_Xhomo','alpha_ak'},true)
    return
end

% Set values from user input, or use defaults
alpha_MK = varg_val(varargin, 'alpha_MK', 95);
alpha_CL = varg_val(varargin, 'alpha_CL', 90);
alpha_Xhomo = varg_val(varargin, 'alpha_Xhomo', 90);
alpha_ak = varg_val(varargin, 'alpha_ak', 95);

% determine if seasons are present
if isstruct(data_tempAgg)
    n=length(data_tempAgg);
    col=fieldnames(data_tempAgg);
    
elseif istimetable(data_tempAgg)
    n=1;
else
    error('data_tempAgg should be either a timetable or a structure of timetable if temporal aggregation is needed (seasonal Mann-Kendall test)');
end

% find the name of the observation in the timetable
if n==1
    obs=fieldnames(data_tempAgg);
elseif n>1
    obs=fieldnames(data_tempAgg(1).(col{1}));
end
% remove the other fields in the timetable:
c=startsWith(obs,["Time"; "Properties";"Variables"]);
obs=obs(~c);
if length(obs)>1
    error('there is more than one observation in the time series');
end
% assemble the temporal aggregations into one time series to compute the
% autocorrelation and the prewhitening
if n==1
    data=data_tempAgg;
    %compute all the prewhitened time series
    % dataPW is a structure containing the 3 prewhitened data.ak_y is the first
    % lag autocorrelation coefficient for the complete time series
    dataPW=prewhite_D(data, (obs{1}), resolution,'alpha_ak',alpha_ak);
    % compute the  Mann-Kendall test without temporal aggregation
    switch PW_method
        case  'PW'
            [result,~ , ~,~ ]=compute_MK_stat(data.Time,dataPW.PW,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
            %%%result.ak=ak_y.PW;
        case  'TFPW_Y'
            [result,~ , ~,~ ]=compute_MK_stat(data.Time,dataPW.TFPW_Y,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
            %%%result.ak=ak_y.TFPW_Y;
        case  'TFPW_WS'
            [result,~ , ~, ~]=compute_MK_stat(data.Time,dataPW.TFPW_WS,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
            %%%result.ak=ak_y.TFPW_WS;
        case  'VCTFPW'
            [result,~ , ~, ~]=compute_MK_stat(data.Time,dataPW.VCTFPW,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
            %%%result.ak=ak_y.VCTFPW;
        case  '3PW'
            [result_PW,~ , ~, ~]=compute_MK_stat(data.Time,dataPW.PW,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
            [result_TFPW_Y,~ , ~, ~]=compute_MK_stat(data.Time,dataPW.TFPW_Y,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
            [result_VCTFPW,~ , ~, ~]=compute_MK_stat(data.Time,dataPW.VCTFPW,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
            %%%result.ak=ak_y.VCTFPW;
            %determine the P and ss
            [result.P, result.ss]=Prob_3PW(result_PW.P,result_TFPW_Y.P, alpha_MK);
            result.slope=result_VCTFPW.slope;
            result.UCL=result_VCTFPW.UCL;
            result.LCL=result_VCTFPW.LCL;
    end
    
elseif n>1
    for m=1:n
        % put the data together to compute the autocorrelation and the prewhitening
        if m==1
            data=data_tempAgg(1).obs;
            for i=2:n
                data=[data;data_tempAgg(i).obs];
            end
            data=sortrows(data);
            %compute all the prewhitened time series
            % dataPW is a structure containing the 3 prewhitened data.ak_y is the first
            % lag autocorrelation coefficient for the complete time series
            dataPW=prewhite_D(data, (obs{1}), resolution,'alpha_ak',alpha_ak);
        end
        
        
        % reassemble the PW time series in temporal aggregation
        ind=data_tempAgg(m).(col{1}).Time;
        dat_mois=dataPW(ind,:);
        
        if sum(~isnan(dat_mois.PW))>1 % un peu trop peu ??
            switch PW_method
                case  'PW'
                    %compute of Mann-Kendall parameters for PW method
                    [result(m), S(m), vari(m), Z(m)]=compute_MK_stat(dat_mois.Time,dat_mois.PW,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
                    %%%ak=ak_y.PW;
                case  'TFPW_Y'
                    %compute of Mann-Kendall parameters for TFPW_Y method
                    [result(m), S(m), vari(m), Z(m)]=compute_MK_stat(dat_mois.Time,dat_mois.TFPW_Y,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
                    %%%ak=ak_y.TFPW_Y;
                case 'TFPW_WS'
                    %compute of Mann-Kendall parameters for TFPW_WS method
                    [result(m), S(m), vari(m), Z(m)]=compute_MK_stat(dat_mois.Time,dat_mois.TFPW_WS,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
                    %%%ak=ak_y.TFPW_WS;
                case  'VCTFPW'
                    %compute of the slope with VCTFPW method
                    [result(m), S(m), vari(m), Z(m)]=compute_MK_stat(dat_mois.Time,dat_mois.VCTFPW,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
                    %%%ak=ak_y.VCTFPW;
                case  '3PW'
                    [result_PW(m),S_PW(m), vari_PW(m), ~]=compute_MK_stat(dat_mois.Time,dat_mois.PW,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
                    [ result_TFPW_Y(m),S_TFPW_Y(m), vari_TFPW_Y(m), ~]=compute_MK_stat(dat_mois.Time,dat_mois.TFPW_Y,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
                    [result_VCTFPW(m), ~, ~, Z_VCTFPW(m)]=compute_MK_stat(dat_mois.Time,dat_mois.VCTFPW,resolution,'alpha_MK',alpha_MK,'alpha_CL',alpha_CL);
                    
                    % write the out file
                    [result(m).P,result(m).ss]=Prob_3PW(result_PW(m).P,result_TFPW_Y(m).P, alpha_MK);
                    
                    % compute the slope and confidence limits for a year [units/y]
                    result(m).slope=result_VCTFPW(m).slope;
                    result(m).UCL=result_VCTFPW(m).UCL;
                    result(m).LCL=result_VCTFPW(m).LCL;
                    % % %                 % compute them as percentage [%/y]
                    % % %                 result(m).median=result_VCTFPW(m).median;
                    % % %                 result(m).slopeP=result_VCTFPW(m).slopeP;
                    % % %                 result(m).UCLP=result_VCTFPW(m).UCLP;
                    % % %                 result(m).LCLP=result_VCTFPW(m).LCLP;
                    % % %                 ak=ak_y.VCTFPW;
                    
            end
            
        else
            
            result(m).P=NaN ;
            result(m).ss= NaN;
            result(m).slope=NaN;
            result(m).UCL=NaN;
            result(m).LCL=NaN;
            % % %         result(m).median=NaN;
            % % %         result(m).slopeP=NaN;
            % % %         result(m).UCLP=NaN;
            % % %         result(m).LCLP=NaN;
        end
    end
    
    %compute for the whole period, that is the whole year
    switch PW_method
        case {'PW','TFPW_Y','TFPW_WS','VCTFPW'}
            %%%result(m+1).ak=ak;
            Ztot=STD_normale_var(nansum(S),nansum(vari));
            if sum(~isnan(dataPW.PW)) > 10
                result(m+1).P=2*(1-normcdf(abs(Ztot),0,1));
            else
                load ('Prob_MK_n');
                result(m+1).P=Prob_MK_n(abs(Stot)+1,sum(~isnan(dataPW.PW)));
            end
            if result(m+1).P<=1-alpha_MK/100
                result(m+1).ss= alpha_MK;
            else
                result(m+1).ss= 0;
            end
            %%%result(m+1).median=nanmedian(dataPW.PW);
            
            %compute xi-carre to test the homogeneity between months.
            %TO DO: change Xhomo as a function of alpha_Xhomo and the number of
            %the seasons
            Xhomo=nansum(Z(1:m).^2)-12*(nanmean(Z(1:m))).^2;
            
        case '3PW'
            % compute the statistical significance for PW
            %%%result(m+1).ak=ak;
            Ztot_PW=STD_normale_var(nansum(S_PW),nansum(vari_PW));
            if sum(~isnan(dataPW.PW)) > 10
                Ptot_PW=2*(1-normcdf(abs(Ztot_PW),0,1));
            else
                load ('Prob_MK_n');
                Ptot_PW=Prob_MK_n(abs(Stot_PW)+1,sum(~isnan(dataPW.PW)));
            end
            % compute the statistical significance for TFPW_Y
            Ztot_TFPW_Y=STD_normale_var(nansum(S_TFPW_Y),nansum(vari_TFPW_Y));
            if sum(~isnan(dataPW.TFPW_Y)) > 10
                Ptot_TFPW_Y=2*(1-normcdf(abs(Ztot_TFPW_Y),0,1));
            else
                load ('Prob_MK_n');
                Ptot_TFPW_Y=Prob_MK_n(abs(Stot_TFPW_Y)+1,sum(~isnan(dataPW.TFPW_Y))); %
            end
            %determine the ss
            [result(m+1).P,result(m+1).ss]=Prob_3PW(Ptot_PW,Ptot_TFPW_Y, alpha_MK);
            
            %compute xi-carre to test the homogeneity between months. Since the slope
            %is computeated from VCTFPW, the homogeneity is also computeated from VCTFPW
            Xhomo=nansum(Z_VCTFPW(1:m).^2)-12*(nanmean(Z_VCTFPW(1:m))).^2;
    end
    %write the yearly slope and CL
    %Xhomo has a chi-squared distributions with n-1 and 1 degree of
    %freedom. Seasonal trends are homogeneous is Xhomo is smaller
    %than the threshold defined by the degree of freedom and the
    %confidence level alpha_Xhomo.
    %change condition: yearly slope not given if the seasons are not
    %homogeneous
    
    if Xhomo<=chi2inv(1-alpha_Xhomo/100,n-1)
        result(m+1).slope=nanmedian([result(1:m).slope]);
        result(m+1).UCL=nanmedian([result(1:m).UCL]);
        result(m+1).LCL=nanmedian([result(1:m).LCL]);
    else
        warning('the trends for the temporal aggregation are not homogeneous');
        result(m+1).slope=NaN;
        result(m+1).UCL=NaN;
        result(m+1).LCL=NaN;
    end
    % %     % result(m+1).median=nanmedian(dataPW.VCTFPW);
    % %     % result(m+1).slopeP=(result(m+1).slope.*100)./abs(result(m+1).median);
    % %     % result(m+1).UCLP=(result(m+1).UCL.*100)'./abs(result(m+1).median);
    % %     % result(m+1).LCLP=(result(m+1).LCL.*100)'./abs(result(m+1).median);
    
end
fclose('all');

