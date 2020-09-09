function [P,ss]=Prob_3PW(P_PW,P_TFPW_Y, alpha_MK)

% estimate the probability of the MK test
% to have the maximal certainty, P is taken as the max of P_PW and P_TFPW_Y
% estimate the statistical significance of the MK test as a function of the
% given confidence level alpha_MK

%in: P_PW: probability computed from the PW prewhitened dataset
%   P_TFPW_Y: probability computed from the TFPW_Y prewhitened dataset
%    alpha_MK: confidence level in % for the MK test. Default=95:

P_alpha=1-alpha_MK/100;
%compute the probability
P=nanmax(P_PW, P_TFPW_Y);
%determine the ss
if P_PW<=P_alpha &&  P_TFPW_Y<=P_alpha
    ss= alpha_MK;
elseif P_PW>P_alpha && P_TFPW_Y<=P_alpha % false positive for TFPW_Y @ alpha %
    ss= -1;
elseif P_TFPW_Y>P_alpha && P_PW<=P_alpha % false positive for TFPW_Y
    ss= -2;
elseif P_TFPW_Y>P_alpha && P_PW>P_alpha % false positive for PW
    ss= 0;
end