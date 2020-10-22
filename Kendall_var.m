% Copyright 2020 MeteoSwiss, contributors of the original matlab version of the code listed in ORIGINAL_AUTHORS
% Copyright 2020 UNIMORE, contributors of the R version of the code listed in AUTHORS
% 
% Distributed under the terms of the BSD 3-Clause License.
% 
% SPDX-License-Identifier: BSD-3-Clause

function y=Kendall_var(data, t,n)
% compute the variance with ties in the data and ties in time
%in:
%       data (array of floats) = data to analyse, must be 1-D
%       t  (array of integer) =  number of element in each tie. Must be 1-D
%       n (array of integer) = nb of non missing data for each year,
%                               correspond to ties in time. Must be 1-D

%out:
%       y (float)= variance

% Source:
%       GAW report No 133 (A. Sirois), p.30 of annexe D

% sanity checks first
if isa(data,'float')==0 || min(size(data))>1
    error('the input "data" of Kendall_var has to be a 1-D array of floats');
end
if isa(t,'float')==0 || min(size(data))>1
    error('the input "t" of Kendall_var has to be a 1-D array of floats');
end
if isa(n,'float')==0 || min(size(data))>1
    error('the input "n" of Kendall_var has to be a 1-D array of floats');
end

% length of data ignoring the NaN
Lreal=sum(~isnan(data));

%compute the variance
var_t1=sum(t.*(t-1).*(2.*t+5));
var_t2=sum(t.*(t-1).*(t-2));
var_t3=sum(t.*(t-1));

var_n1=sum(n.*(n-1).*(2.*n+5));
var_n2=sum(n.*(n-1).*(n-2));
var_n3=sum(n.*(n-1));

L=max(size(data));
y=(((Lreal *(Lreal-1)*(2*Lreal+5))-var_t1-var_n1)/18)+(var_t2 *var_n2/(9*Lreal*(Lreal-1)*(Lreal-2)))+(var_t3 *var_n3/(2*Lreal*(Lreal-1)));
fclose('all');