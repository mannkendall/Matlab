% varg_proof v1.1: licensed by Tami Bond yark@colostate.edu under Creative Commons by-sa 2.0

%% varg_proof(vargin, oknames, giveoutput) is an operator for varargin. It checks
%%   whether all the arguments are present in oknames.

function isok = varg_proof(vargin, oknames, giveoutput)

   if nargin<3
       giveoutput = false;
   end

   badnames = '';
   varnames = vargin(1:2:length(vargin));

   % check that all varnames are in oknames
   for i=1:length(varnames)
       if sum(strcmp(varnames(i), oknames))==0
           badnames = [badnames cell2mat(varnames(i)) ' '];
       end
   end
   
   isok = (length(badnames)==0);
   if ~isok & giveoutput
       disp(['Invalid arguments:' badnames]);
   end
       
return
   
  