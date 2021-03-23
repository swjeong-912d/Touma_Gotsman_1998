function [vertex,face] = read_off(filename)

% read_off - read data from OFF file.
%
%   [vertex,face] = read_off(filename);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) 2003 Gabriel Peyr?


fid = fopen(filename,'r');
if( fid==-1 )
    warning(['Can''t open the file: ', filename]);
    vertex = zeros(100000,1);
    face = zeros(100000,1);
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
[a,str] = strtok(str); nvert = str2num(a);
[a,str] = strtok(str); nface = str2num(a);
A = textscan(fid,'%f %f %f', nvert);
if ~isequal([length(A) length(A{1})], [3 nvert])
    warning('Problem in reading vertices.');
end
vertex = cell2mat(A);
A = textscan(fid,'%d %d %d %d', nface);
if ~isequal([length(A) length(A{1})], [4 nface])
    warning('Problem in reading faces.');
end
face = double(cell2mat(A(:,2:4))+1);


fclose(fid);

