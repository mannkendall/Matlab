%  Copyright 2020 MeteoSwiss, contributors of the original matlab version of the code listed in ORIGINAL_AUTHORS
% 
% Distributed under the terms of the BSD 3-Clause License.
% 
% SPDX-License-Identifier: BSD-3-Clause

function t=Nb_tie(data,resolution)
 
 % compute the number of data considered equivalent and treated as ties
 
 %in:
 %      data (array of floats)= data to be analysed, must be 1-D
 %      resolution ((float)= the measurement resolution, i.e. delta value
 %      below which 2 measurements are considered equivalent
 
 %out: 
 %      t (n arrays of int)= amount of ties in the data 
 
 % check if data is a one dimentional array of floatss
 if isa(data,'float')==0 || min(size(data))>1
     error('the input "data" of Nb_tie_D has to be a 1-D array of floats');
 end
 %check if resolution if a float
  if isa(resolution,'float')==0 || max(size(resolution))>1
     error('the input "resolution" of Nb_tie_D has to be a single float value');
 end
 
 if sum(~isnan(data))>4
     m=min(data);
     M=max(data);
     %determine the bin edges containing the data
     if (m<0 & M<0) | (m>0 & M>0)
         interval=abs(M-m);
     else
         interval=M+abs(m);
     end
     if  resolution>=interval
         error('the given resolution is too large for the considered dataset');
     end
     step=floor(interval/resolution)+1;
     % if the time series is too big and generate a memory error
     % the number of bins can be decreased with usually low impact on the Mann-Kendall results
% %      while step>1000
% %          resolution=resolution*2;
% %          step=interval/resolution;
% %      end
     t=histcounts(real(data),'BinLimits',[nanmin(data),nanmin(data)+step*resolution],'BinWidth',resolution)'; 
 else %if the time series contains less than 4 data, return NaN
     t=NaN;
 end
fclose('all');
