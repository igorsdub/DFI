function [ ] = testfunc(list)
%UNTITLED Test Octave 
%   Detailed explanation goes here
hlist = list; 
l1 = [1,11,22,33,44,55];
l2 = l1(hlist);
outfile = fopen('test.dat','w');
for i=l2
    fprintf(outfile,'%d\n',i);
end
fclose(outfile);
end

