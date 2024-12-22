function [x,y]=meshline(lnr,dofn)
% assign nodnumbers for a line

global CALFEM_COORD
global CALFEM_DOF
global CALFEM_LINES
global CALFEM_LINES_NODES
global CALFEM_POINTS
global CALFEM_POINTS_NODE
if size(CALFEM_LINES_NODES,1)>0     % check if done before
   if find(CALFEM_LINES_NODES(:,1)==lnr)
      [x,y]=divline(lnr);
      break
   end
end
CALFEM_LINES_NODES(lnr,1)=lnr;
nextnode=size(CALFEM_COORD,1)+1;
[nodes,dfsb]=size(CALFEM_DOF);
if dofn>dfsb     % more dofs per node than before
   CALFEM_DOF=[CALFEM_DOF zeros(nodes,dofn-dfsb)];
end
if size(CALFEM_POINTS_NODE,1)>0     % nodes done before  
   indx=find(CALFEM_POINTS_NODE(:,1)==CALFEM_LINES(lnr,3));
   if indx       % first node already been assigned
      CALFEM_LINES_NODES(lnr,2)=CALFEM_POINTS_NODE(indx,2);
   else
      CALFEM_POINTS_NODE(CALFEM_LINES(lnr,3),:)=[CALFEM_LINES(lnr,3) nextnode];
      CALFEM_LINES_NODES(lnr,2)=nextnode;
      CALFEM_COORD=[CALFEM_COORD; CALFEM_POINTS(CALFEM_LINES(lnr,3),2:4)];
      CALFEM_DOF=[CALFEM_DOF; (nextnode-1)*dofn+1:nextnode*dofn...
            zeros(1,dfsb-dofn)];
      nextnode=nextnode+1;
   end
   % the nodes in between
   [x,y]=divline(lnr);
   for k=2:(CALFEM_LINES(lnr,2))
      CALFEM_LINES_NODES(lnr,k+1)=nextnode;
      CALFEM_COORD=[CALFEM_COORD; x(k) y(k) 0];
      CALFEM_DOF=[CALFEM_DOF; (nextnode-1)*dofn+1:nextnode*dofn...
            zeros(1,dfsb-dofn)];
      nextnode=nextnode+1;
   end
   % last node
   indx=find(CALFEM_POINTS_NODE(:,1)==CALFEM_LINES(lnr,4));
   if indx      % last node already been assigned
      CALFEM_LINES_NODES(lnr,CALFEM_LINES(lnr,2)+2)=CALFEM_POINTS_NODE(indx,2);
   else
      CALFEM_POINTS_NODE(CALFEM_LINES(lnr,4),:)=[CALFEM_LINES(lnr,4) nextnode];
      CALFEM_LINES_NODES(lnr,CALFEM_LINES(lnr,2)+2)=nextnode;
      CALFEM_COORD=[CALFEM_COORD; CALFEM_POINTS(CALFEM_LINES(lnr,4),2:4)];
      CALFEM_DOF=[CALFEM_DOF; (nextnode-1)*dofn+1:nextnode*dofn...
            zeros(1,dfsb-dofn)];
   end
else
   [x,y]=divline(lnr);
   for k=1:(CALFEM_LINES(lnr,2)+1)
      CALFEM_LINES_NODES(lnr,k+1)=k;
      CALFEM_COORD=[CALFEM_COORD; x(k) y(k) 0];
      CALFEM_DOF=[CALFEM_DOF; (k-1)*dofn+1:k*dofn];
   end
   CALFEM_POINTS_NODE(CALFEM_LINES(lnr,3),:)=[CALFEM_LINES(lnr,3) 1];
   CALFEM_POINTS_NODE(CALFEM_LINES(lnr,4),:)=[CALFEM_LINES(lnr,4) k];
end      
 
