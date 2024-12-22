function gnew(typ)
% gnew(typ)
%------------------------------------------------------------------------------
% PURPOSE
%  Delete part of the geometry
%
% INPUT: typ  type of geometry to be deleted ('point','line',surf' 
%             ,'mesh' or 'geom')
%
%-------------------------------------------------------------

% REFERENCES
%   M Ristinmaa 1994-05-24
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology

global CALFEM_POINTS
global CALFEM_LINES
global CALFEM_SURFS
global CALFEM_POINTS_NODE
global CALFEM_LINES_NODES

typs=typ(1:3);
dm=0;
if typs=='poi' | typs=='geo'
  clf
  disp('All points deleted ')
  disp('All lines deleted ')
  disp('All surfaces deleted ')
  clear global CALFEM_POINTS
  clear global CALFEM_LINES
  clear global CALFEM_SURFS
  clear global CALFEM_POINTS_NODE
  clear global CALFEM_LINES_NODES
  dm=1;
elseif typs=='lin'
  clf
  disp('All lines deleted ')
  disp('All surfaces deleted ')
  clear global CALFEM_LINES
  clear global CALFEM_SURFS
  clear global CALFEM_POINTS_NODE
  clear global CALFEM_LINES_NODES
  gdraw('poi')
  dm=1;
elseif typs=='sur'
  clf
  disp('All surfaces deleted ')
  clear global CALFEM_SURFS
  clear global CALFEM_POINTS_NODE
  clear global CALFEM_LINES_NODES
  gdraw('poi'); gdraw('line');
  dm=1;
end  
if typs=='mes' | dm
  disp('Mesh deleted')
  clear global CALFEM_COORD
  clear global CALFEM_EDOF
  clear global CALFEM_DOF
  clear global CALFEM_LINES_NODES
  clear global CALFEM_POINTS_NODE
  clear global CALFEM_LINES_ELEMS
  clear global CALFEM_SURFS_ELEMS
  clf
  gdraw('poi'); gdraw('line');
end
