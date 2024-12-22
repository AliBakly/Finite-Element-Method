function [f]=gcload(typ,nr,dnr,val)
% [f]=gcload(typ,nr,dnr,val)
%------------------------------------------------------------------------------
% PURPOSE
%  Apply concentrated loads on the defined geometry
%
% INPUT:  typ  defines the part to which 'force is mapplied ('point' or 'line')
%
%         nr   number of the part
%
%         dnr  the dof to which the force is applied
%
%         val  value of the force
%-------------------------------------------------------------

% REFERENCES
%   M Ristinmaa 1994-05-24
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
global CALFEM_POINTS_NODE
global CALFEM_LINES_NODES
global CALFEM_LINES
global CALFEM_DOF
tp=typ(1:3);
dofs=max(max(CALFEM_DOF));
f=zeros(dofs,1);
if tp=='poi'
   pl=find(CALFEM_POINTS_NODE(:,1)==nr);
   if isempty(pl)
      disp('GEOMETRY ERROR: point not yet assigned dofs')
      return;
   end
   f(CALFEM_DOF(CALFEM_POINTS_NODE(pl,2),dnr))=val;
elseif tp=='lin'
   ll=find(CALFEM_LINES(:,1)==nr);
   if CALFEM_LINES(ll,2)~=0   % not a line combine
      if CALFEM_LINES_NODES(ll,1)==0
         disp('GEOMETRY ERROR: line not yet assigned dofs')
         return;
      end
      nds=max(find(CALFEM_LINES_NODES(ll,:)>0));
      for k=2:nds
         f(CALFEM_DOF(CALFEM_LINES_NODES(ll,k),dnr))=val;
      end
   else
      f=cload('lin',CALFEM_LINES(nr,3),dnr,val);
      f1=cload('lin',CALFEM_LINES(nr,4),dnr,val);
      f=f+f1;
   end
end
