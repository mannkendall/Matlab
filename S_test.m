% Copyright 2020 MeteoSwiss, contributors of the original matlab version of the code listed in ORIGINAL_AUTHORS
% 
% Distributed under the terms of the BSD 3-Clause License.
% 
% SPDX-License-Identifier: BSD-3-Clause

function [S, n]=S_test(data, time)
% compute the S statistic (Si) for the MK test

%input:    
%       data (array of floats)= data. Must be 1-D
%       time (array of float = time. Must be given in datenum(time)

%out    S (float)= double sum on the sign of the difference between data pairs (Si)      
%       n (array of integer)= number of valid data in each year of the time series

%source:
%       Gilbert, 1987

% sanity checks first
if isa(data,'float')==0 || min(size(data))>1
    error('the input "data" of S_test has to be a 1-D array of floats');
end
if isa(time,'float')==0 || min(size(time))>1
    error('the input "time" of S_test has to be a 1-D array of floats');
end
if length(data)~=length(time)
    error(' data and time should have the same length');
end
% put the time in a vector
t=datevec(time);
an_debut=nanmin(t(:,1));
an_fin=nanmax(t(:,1));
nb_an=an_fin-an_debut+1;
% initialise Sij
Sij=NaN(nb_an,nb_an);

for i=an_debut:an_fin
    %compute the number of valid data in each year
    n(i-an_debut+1,1)=sum(~isnan(data(t(:,1)==i,1)));
    %compute the number of positive and negative pair differences
    if i<an_fin
        for j=i+1:an_fin
            nj=size(data(t(:,1)==j),1);
            sj=ones(nj,1);
            sj(:)=NaN;
            
            dataj=data(t(:,1)==j,1);
            for k=1:nj
                sj(k)=nansum(sign(dataj(k)-data(t(:,1)==i,1)));
            end
            Sij(j-an_debut+1,i-an_debut+1)=nansum(sj);
        end
    end
end 
S=nansum(nansum(Sij));
fclose('all');
