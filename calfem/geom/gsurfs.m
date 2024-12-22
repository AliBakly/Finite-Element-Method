function s=gsurfs(snum,l1,l2,l3,l4)
% s=gsurfs(snum,l1,l2,l3,l4)
%------------------------------------------------------------------------------
% PURPOSE
%  Define a surface which is used to describe the geometry
%
% INPUT:  snum        number of surface, if 0 a number will be generated
%
%         l1,l2,l3,l4 border lines of surface
%
% OUTPUT: s           if no arguments the list of surfaces will be returned
%-------------------------------------------------------------

% REFERENCES
%   M Ristinmaa 1994-05-24
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology

global CALFEM_SURFS;
if nargin==0
   s=CALFEM_SURFS;
   return
end
if nargin==4 l4=0; end
nrsurfs=size(l1,1);
if nrsurfs>1     % many surfaces specified in a matrix
   for i=1:nrsurfs
      if size(snum,1)<i snum(i)=0; end
      if size(l4,1)<i l4(i)=0; end
      gsurfs(snum(i),l1(i),l2(i),l3(i),l4(i));
   end
   break
end

[r,c]=size(CALFEM_SURFS);
if r==0
   if snum==0
      snum=1;
   end
   CALFEM_SURFS(snum,1:5)=[snum l1 l2 l3 l4];
else
   if snum==0
      n=max(CALFEM_SURFS(:,1));
      snum=n+1;
   end
   if isempty(find(CALFEM_SURFS(:,1)==snum))
      CALFEM_SURFS(snum,:)=[snum l1 l2 l3 l4];
   else
      disp('GEOMETRY ERROR: surfacenumber already in use')
   end
end
