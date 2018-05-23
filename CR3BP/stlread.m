function geometry = stlread(filename)

%STLREAD reads an stl file in binary format into vertex and face matrices.
%   [V,F,N,C,STLTITLE] = STLREAD(FILENAME) reads an stl file in binary 
%   format with specified FILENAME and gives in output matrices containing
%   the vertices for all triangles V, vertex lists defining each triangle
%   face F, normals for each triangle face N and 5 bits RGB color data for
%   each triangle face C. STLTITLE contains the title of the specified stl
%   file.
%
%   To see plot the 3D surface use:
%       patch('Faces',F,'Vertices',V,'FaceVertexCData',C);
%   or
%       plot3(V(:,1),V(:,2),V(:,3),'.');
%
%   Duplicate vertices can be removed using:
%       [V,F] = patchslim(V,F);
%
%   For more information see:
%   http://www.esmonde-white.com/home/diversions/matlab-program-for-loading-stl-files
%
%   Based on code originally written by:
%       Doron Harlev
%   and combined with some code by:
%       Eric C. Johnson, 11-Dec-2008
%   Copyright 1999-2008 The MathWorks, Inc.
% 
%   Re-written and optimized by Francis Esmonde-White, May 2010.

use_color = nargout>=4;

fid = fopen(filename,'r'); %Open the file, assumes STL Binary format.
if fid == -1 
    error('File could not be opened, check name or path.')
end

ftitle = fread(fid,80,'uchar=>schar'); % Read file title
numFaces = fread(fid,1,'int32'); % Read number of Faces

T = fread(fid,inf,'uint8=>uint8'); % read the remaining values
fclose(fid);

stltitle = char(ftitle');

% Each facet is 50 bytes
%  - Three single precision values specifying the face normal vector
%  - Three single precision values specifying the first vertex (XYZ)
%  - Three single precision values specifying the second vertex (XYZ)
%  - Three single precision values specifying the third vertex (XYZ)
%  - Two color bytes (possibly zeroed)

% 3 dimensions x 4 bytes x 4 vertices = 48 bytes for triangle vertices
% 2 bytes = color (if color is specified)

trilist = 1:48;

ind = reshape(repmat(50*(0:(numFaces-1)),[48,1]),[1,48*numFaces]) + repmat(trilist,[1,numFaces]);
tri = reshape(typecast(T(ind),'single'),[3,4,numFaces]);

n = squeeze(tri(:,1,:))';
n = double(n);

v = tri(:,2:4,:);
v = reshape(v,[3,3*numFaces]);
v = double(v)';

f = reshape(1:3*numFaces,[3,numFaces])';

c0 = typecast(T(49:50),'uint16');
if bitget(c0(1),16) == 1
    trilist = 49:50;
    ind = reshape(repmat(50*(0:(numFaces-1)),[2,1]),[1,2*numFaces])+repmat(trilist,[1,numFaces]);
    c0 = reshape(typecast(T(ind),'uint16'),[1,numFaces]);
    r = bitshift(bitand(2^16-1,c0),-10);
    g = bitshift(bitand(2^11-1,c0),-5);
    b = bitand(2^6-1,c0);
    c = [r;g;b]';
else
    c = zeros(numFaces,3);
end

geometry = struct('Vertices',v,'Faces',f,'Normals',n,'FaceColor',c,'Title',stltitle);

end