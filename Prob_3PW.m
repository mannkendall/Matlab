function [P,ss]=Prob_3PW(P_PW,P_TFPW_Y, alpha_MK)

% 1) Estimate the probability of the MK test with the 3PW method.
%    To have the maximal certainty, P is taken as the maximum of P_PW and P_TFPW_Y
% 2) Estimate the statistical significance of the MK test as a function of the
%    given confidence level alpha_MK

%input: 
%       P_PW (float): probability computed from the PW prewhitened dataset
%       P_TFPW_Y (float): probability computed from the TFPW_Y prewhitened dataset
%       alpha_MK (float): confidence level in % for the MK test. 

% Some sanity checks first
if isa(P_PW,'float')==0 || min(size(P_PW))>1
    error('the input "P_PW" of Prob_3PW has to be a float');
end
if isa(P_TFPW_Y,'float')==0 || min(size(P_TFPW_Y))>1
    error('the input "P_TFPW_Y" of Prob_3PW has to be a float');
end
if isa(alpha_MK,'float')==0 || min(size(alpha_MK))>1
    error('the input "alpha_MK" of Prob_3PW has to be a float');
end

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