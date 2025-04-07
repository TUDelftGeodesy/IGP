function stmprintopt(value,indent)
%stmprintopt  Pretty print STM options. 
%  STMPRINTOPT(VALUE) prints whatever is in VALUE (char, string, array, 
%  cell, struct) recursively in a nice manner.
% 
%  The function can be called recursively with a second INDENT argument.
%
%  (c) Hans van der Marel, Delft University of Technology, 2022.

% Created:  26 October 2022 by Hans van der Marel
% Modified: 7 March 2024 by Hans van der Marel
%             - added missing quotes to cellstr output after ;
%             - fix output of structure arrays

if nargin < 2
    indent='';
end

if isempty(value)
   fprintf('<empty>')
elseif ischar(value) || isstring(value) || isnumeric(value) || islogical(value)
   if numel(value) < 1
      fprintf('<empty>')
   else
      fprintf('%s', mat2str(value))
   end
elseif iscellstr(value)
   tmp='{ ''';
   for m=1:size(value,1)
      tmp = [ tmp strjoin(value(m,:),"' '") ''' ; '''];
   end
   tmp(end)='}';
   fprintf('%s',tmp) 
elseif iscell(value)
   if numel(value) < 1
      fprintf(' <empty>')
   else
      % call stmprintopt recursively
      fprintf('{');
      for m=1:numel(value)
          stmprintopt(value{m})
      end
      fprintf('}')
   end
elseif isstruct(value)
   if numel(value) < 1
      fprintf(' <empty>')
   else
      for m=1:numel(value)
          fields=fieldnames(value(m));
          if numel(fields) < 1
             fprintf(' <empty>')
          else
             % call stmprintopt recursively
             for k=1:numel(fields)
                if k > 1 || numel(indent) > 1 fprintf('\n'); end
                fprintf('%s%s: ',indent,fields{k})
                stmprintopt(value(m).(fields{k}),[indent '   '])
             end
          end
      end
   end
else
   disp(value)
end

end