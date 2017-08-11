function aus = getfilenamenumber(filename)

%prefix_length     = strfind(setup.akt_filename,'_') ; if isempty(prefix_length) prefix_length=length(setup.akt_filename) ; else prefix_length = max([ (prefix_length-1)  13 ]) ; end ;
pos                = find(filename<=57&filename>=48) ;      % ASCII  (0-9) -> Decimal (48-57)
help               = pos(1) ;
while any(help==pos)  
    help = help+1 ;  
end ;
pos                = (pos(1):help-1) ;

aus = filename(pos) ;