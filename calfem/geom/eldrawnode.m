function eldrawnode(edof,coord,Dof,nen,fromnode)
% eldrawnode(edof,coord,Dof,nen)
% draw mesh with node numbers

if nargin<5
   fromnode=1;
end
[ex,ey]=coordxtr(edof,coord,Dof,nen);
eldraw2(ex,ey,[1 4 0]);
[nnodes dpn]=size(Dof);
delta=max(max(ex))/100;
for i=fromnode:nnodes
   text(coord(i,1)+delta,coord(i,2)+delta/2,int2str(i));
end