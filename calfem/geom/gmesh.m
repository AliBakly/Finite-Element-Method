function [edof,coord,Dof]=gmesh(type,nbr,eltype,ndn,drawop)
% [edof,coord,Dof]=gmesh(type,nr,eltyp,ndn,drawop)
%----------------------------------------------------------------
%  PURPOSE
%   Generates elements on a given geometry
%
%  INPUT:  type    describes the part to be meshed ('line' or 'surf')
%
%          nr      is number of the part to be meshed
%
%          eltype  is the element type which will be used in mesh
%                    'be2' two nodes beam element
%                    'tr3' three node triangular surface element
%                    'tr6' six node triangular surface element
%                    'qu4' four   "      "      "
%                    'qu8' eight  "      "      "
%                    'qr3' three node triangular surface element in rect. conf.
%                    'qr6' six node triangular surface element in rect. conf.
%
%          ndn     is number of dofs per node
%
%          drawop  controls drawing of mesh
%                    0    don't draw mesh
%                    1    draw mesh without elementnumber (default)
%                    2      "    "  with        "
%		     3      "    "    "   nodenumbers
%		     4      "    "    "        "      and elementnumbers
%
%  OUTPUT: edof    topology matrix for meshed areas
%
%          coord   global coordinate matrix 
%
%          Dof     global nodal dof matrix
%
% -------------------------------------------------------------

% REFERENCES
%   A Magnusson 1995-09-07
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology

global CALFEM_POINTS
global CALFEM_LINES
global CALFEM_SURFS
global CALFEM_COORD
global CALFEM_EDOF
global CALFEM_DOF
global CALFEM_POINTS_NODE
global CALFEM_LINES_NODES
global CALFEM_LINES_ELEMS
global CALFEM_SURFS_ELEMS

if nargin==4
   drawop=1;   %draw mesh
end

nodes=size(CALFEM_DOF,1);
if type(1)=='l'
   if eltype=='be2'   	                  npel=2;
   elseif eltype=='be3'			  npel=3;
   else 
      disp('GEOMETRY ERROR: element type not implemented');
      err=1;
      break
   end
   err=linemesh(nbr,eltype,ndn,npel);
elseif type(1)=='s'
   if eltype=='tr3'|eltype=='qr3'   	  npel=3;
   elseif eltype=='qu4'			  npel=4;
   elseif eltype=='qu8'			  npel=8;
   elseif eltype=='qr6'			  npel=6;
   else 
      disp(['Sorry, element ',eltype,' has not yet been implemented'])
      err=1;
      break
   end
   if (npel==8 | npel==6) & ...
      (rem(CALFEM_LINES(CALFEM_SURFS(nbr,2),2),2) | rem(CALFEM_LINES(CALFEM_SURFS(nbr,3),2),2))
      disp(['For element type ',eltype,' the lines have to be devided int an even number']);
      err=1;
      break
   end
   err=meshsurf(nbr,eltype,ndn,npel);
end
if drawop>0 & ~err
   [edof,coord,Dof]=getgeom(type,nbr);
   [ex,ey]=coordxtr(edof,coord,Dof,npel);
   if drawop==1
      eldraw2(ex,ey,[1 5 3]);
   elseif drawop==2
      eldraw2(ex,ey,[1 5 3],edof(:,1));
   elseif drawop==3
      eldrawnode(edof,coord,Dof,npel,nodes+1);
   elseif drawop==4
      eldraw2(ex,ey,[1 5 3],edof(:,1));
      eldrawnode(edof,coord,Dof,npel,nodes+1);
   end
end
if nargout>0 & ~err
   [edof,coord,Dof]=getgeom(type,nbr);
else
   edof=[];
end
