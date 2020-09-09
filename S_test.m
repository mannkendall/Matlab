function [S, n]=S_test(data, time)
% compute the S statistic (Si) for the MK test

%in:    data= contains data to be analyzed
%       time = datenum(time)

%out    S is the double sum on the sign of the difference between data (Si)      
%       n=number of valid data in each year of the time series

if length(data)~=length(time)
    error(' data and time should have the same length');
end

t=datevec(time);
an_debut=nanmin(t(:,1));
an_fin=nanmax(t(:,1));
nb_an=an_fin-an_debut+1;
colonne=1;
Sij=NaN(nb_an,nb_an);

for i=an_debut:an_fin
    n(i-an_debut+1,1)=sum(~isnan(data(t(:,1)==i,1)));
    if i<an_fin
        for j=i+1:an_fin
            nj=size(data(t(:,1)==j),1);
            sj=ones(nj,1);
            sj(:)=NaN;
            
            dataj=data(t(:,1)==j,colonne);
            for k=1:nj
                sj(k)=nansum(sign(dataj(k)-data(t(:,1)==i,colonne)));
            end
            Sij(j-an_debut+1,i-an_debut+1)=nansum(sj);
        end
    end
end 
S=nansum(nansum(Sij));
fclose('all');
