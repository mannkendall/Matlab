function y=Kendall_var(data, t,n)
% calculate the variance with ties in the data and ties in time 
%in:    data = one array of data
%       t = the number of element in each tie 
%       n= nb of non missing data for each year of the time series

%out: y= variance
% taken from GAW report No 133 (A. Sirois), p.30 of annexe D

%%% put a test to verify that data, t and n have the right format
%%% put a test to verify that t and n have the same size

Lreal=sum(~isnan(data));

var_t1=sum(t.*(t-1).*(2.*t+5));
var_t2=sum(t.*(t-1).*(t-2));
var_t3=sum(t.*(t-1));

var_n1=sum(n.*(n-1).*(2.*n+5));
var_n2=sum(n.*(n-1).*(n-2));
var_n3=sum(n.*(n-1));

L=max(size(data));
y=(((Lreal *(Lreal-1)*(2*Lreal+5))-var_t1-var_n1)/18)+(var_t2 *var_n2/(9*Lreal*(Lreal-1)*(Lreal-2)))+(var_t3 *var_n3/(2*Lreal*(Lreal-1)));
 fclose('all');