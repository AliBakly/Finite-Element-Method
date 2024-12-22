function [edof,coord,Dof]=getgeom(typ,nr)
% [edof,coord,Dof]=getgeom(typ,nr)
%------------------------------------------------------------------------------
% PURPOSE
%  Extract the topology and coordinates from the defined geometry
%
% INPUT:  typ    describes what kind of parts geometry to get ('line' or 'surf')
%
%         nr     number of part 
%
% OUTPUT: coord  global coordinate matrix 
%
%         edof   topology martix
%
%         Dof    global nodal dof matrix
%
%-------------------------------------------------------------

% REFERENCES
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology

global CALFEM_COORD
global CALFEM_EDOF
global CALFEM_DOF
global CALFEM_LINES_ELEMS
global CALFEM_SURFS_ELEMS
coord=CALFEM_COORD;
Dof=CALFEM_DOF;
if nargin==0
   edof=CALFEM_EDOF;
   break
end
tp=typ(1:3);
if tp=='lin'
   lnr=find(CALFEM_LINES_ELEMS(:,1)==nr);
   fr=find(CALFEM_EDOF(:,1)==CALFEM_LINES_ELEMS(lnr,2));
   tr=find(CALFEM_EDOF(:,1)==CALFEM_LINES_ELEMS(lnr,3));
   mx=max(find(CALFEM_EDOF(fr,:)>0));
   edof=CALFEM_EDOF(fr:tr,1:mx);
elseif tp=='sur'
   lnr=find(CALFEM_SURFS_ELEMS(:,1)==nr);
   fr=find(CALFEM_EDOF(:,1)==CALFEM_SURFS_ELEMS(lnr,2));
   tr=find(CALFEM_EDOF(:,1)==CALFEM_SURFS_ELEMS(lnr,3));
   mx=max(find(CALFEM_EDOF(fr,:)>0));
   edof=CALFEM_EDOF(fr:tr,1:mx);
end
