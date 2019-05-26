function out = print_doc(name_doc)
% PRINT_DOC  prints in a pretty way the docs for a field
%   Prints argument in specific way for docs.
global s

try
  description = eval(name_doc);
  description = strrep(description,'\n','\n\t');
catch exception
  description = 'needs documentation';
end
name_doc = regexprep(name_doc,'..docs\.','.');
disp(sprintf([name_doc,'\n\t',description]))
end