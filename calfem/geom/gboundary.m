function [bc]=gboundary(typ,nr,dnr,val,bc)
% [bc]=gboundary(typ,nr,dnr,val,bc)
%----------------------------------------------------------------
% PURPOSE
%   Puts boundary condition onto the defined geometry
%
% INPUT:  typ  defines the part to which boundary condition is
%              applied ('point' or 'line')
%
%         nr   number of the part
%
%         dnr  the dof to which the bc is applied
%
%         val  value to be prescribed (default is zero)
% 
%         bc   if bc is given as input, it will be updated
%
% OUTPUT: bc   the boundary condition matrix
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
if nargin==3 val=0; bc=[]; end
if nargin==4 bc=[]; end
tp=typ(1:3);
done=0;
if tp=='poi'
   pl=find(CALFEM_POINTS_NODE(:,1)==nr);
   if isempty(pl)
      disp('GEOMETRY ERROR: point not yet assigned dofs') 
      return;
   end
   if ~isempty(bc)
      done=0;
      if find(bc(:,1)==CALFEM_DOF(CALFEM_POINTS_NODE(pl,2),dnr))
         done=1;
      end
   end
   if ~done
      bc=[bc; CALFEM_DOF(CALFEM_POINTS_NODE(pl,2),dnr) val];
   end
elseif tp=='lin'
   ll=find(CALFEM_LINES(:,1)==nr);
   if CALFEM_LINES(ll,2)~=0   % not a line combine
      if CALFEM_LINES_NODES(ll,1)==0
         disp('GEOMETRY ERROR: line not yet assigned dofs') 
         return;
      end
      nds=max(find(CALFEM_LINES_NODES(ll,:)>0));
      for k=2:nds
         done=0;
         if ~isempty(bc)
            if find(bc(:,1)==CALFEM_DOF(CALFEM_LINES_NODES(ll,k),dnr))
               done=1;
            end
         end
	 if ~done
            bc=[bc; CALFEM_DOF(CALFEM_LINES_NODES(ll,k),dnr) val];
	 end
      end
   else
      bc=gboundary('lin',CALFEM_LINES(nr,3),dnr,val,bc);
      bc=gboundary('lin',CALFEM_LINES(nr,4),dnr,val,bc);
   end
end
