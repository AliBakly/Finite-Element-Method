function p=gpoints(pnum,x,y,z,dr)
% p=gpoints(pnum,x,y,z,dr)
%------------------------------------------------------------------------------
% PURPOSE
%  Create a point which is used to define the geometry
%  
% INPUT:  pnum     number of point, if 0 a number will be generated
%                  If first argument is a matrix points will create
%                  all the points specified (as long as no nubers clash)
%
%         x, y, z  coordinates of point (y and z default =0)
%
%         dr       if the point is not to be drawn dr=0
%
% OUTPUT: p        if no arguments the list of points will be returned
%-------------------------------------------------------------

% REFERENCES
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
global CALFEM_POINTS;
if nargin==0
   p=CALFEM_POINTS;
   return
end
if nargin==4 dr=1; end
if nargin==3 z=0; dr=1; end
if nargin==2 y=0; z=0; dr=1; end
nrpoints=max(size(pnum,1),size(x,1));
if nrpoints>1       % many points specified in a matrix
   for i=1:nrpoints
      if size(pnum,1)<i pnum(i)=0; end
      if size(x,1)<i x(i)=0; end
      if size(y,1)<i y(i)=0; end
      if size(z,1)<i z(i)=0; end
      if size(dr,1)<i dr(i)=1; end
      gpoints(pnum(i),x(i),y(i),z(i),dr(i));
   end
   break
end
[r,c]=size(CALFEM_POINTS);
if r==0
   if pnum==0
      pnum=1;
   end
   CALFEM_POINTS(pnum,1:4)=[pnum x y z];
   if dr gdraw('poi',1,pnum); end
else
   if pnum==0
      n=max(CALFEM_POINTS(:,1));
      pnum=n+1;
      disp(['Pointnumber ',num2str(pnum),' created'])
   end
   if isempty(find(CALFEM_POINTS(:,1)==pnum))
      CALFEM_POINTS(pnum,:)=[pnum x y z];
      if dr gdraw('poi',1,pnum); end
   else
      disp('GEOMETRY ERROR: pointnumber already in use')
   end
end
