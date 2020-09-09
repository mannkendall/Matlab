 function t=Nb_tie_D(data,resolution)
 
 % compute the number of data considered as equivalent and treated as ties
 
 %in: data= data to be analysed

 % resolution is the measurement resolution, i.e. the interval between which 2 measures are considered as equals
 
 % out: vector with the number of data per tie 
 
 if sum(~isnan(data))>4
     m=min(data);
     M=max(data);
     %determine the domains containing the data
     if (m<0 & M<0) | (m>0 & M>0)
         interval=abs(M-m);
     else
         interval=M+abs(m);
     end
     if  resolution>=interval
         error('the given resolution is too large for the considered dataset');
     end
     step=interval/resolution;
     while step>1000
         resolution=resolution*2;
         step=interval/resolution;
     end
     t=hist(real(data),step)';
     fclose('all');
 else
     t=NaN;
 end

