function [X,Y,nodes]=surfcoords(borders,eltype,dofn)
% calculate coordinates for all nodes in surface

global CALFEM_COORD
global CALFEM_LINES_NODES
global CALFEM_LINES
global CALFEM_DOF

[x1,y1]=meshline(borders(1),dofn);
[x2,y2]=meshline(borders(2),dofn);
[x3,y3]=divline(borders(3)); % don't mesh this line yet
[x4,y4]=meshline(borders(4),dofn);


ny=size(x2,2); nx=size(x1,2);
nextnode=size(CALFEM_COORD,1)+1;
dfsb=size(CALFEM_DOF,2);

X=zeros(ny,nx); Y=X; nodes=X;
% coordinates on border
X(1,:)  =x1; Y(1,:)  =y1;                              % line 1
X(ny,:) =x3; Y(ny,:) =y3;			       % line 3
X(2:ny-1,nx) =x2(2:ny-1)'; Y(2:ny-1,nx) =y2(2:ny-1)';  % line 2
X(2:ny-1,1)  =x4(2:ny-1)'; Y(2:ny-1,1)  =y4(2:ny-1)';  % line 4
% nodenumbers on border
nodes(1,:)    =CALFEM_LINES_NODES(borders(1),2:nx+1);
nodes(1:ny,nx)=CALFEM_LINES_NODES(borders(2),2:ny+1)';
nodes(1:ny,1) =CALFEM_LINES_NODES(borders(4),2:ny+1)';

% check for curved lines
culs=0; cul=zeros(4,1);
for k=1:4
   if CALFEM_LINES(borders(k),5)~=0
      culs=culs+1;
      cul(k)=1;
   end
end
if culs
   [X,Y,nodes]=curvedlines(cul,borders,X,Y,nodes,dofn);
else
   % internal nodes
   for k=2:nx-1       % kolonn
      for j=2:ny-1    % row
       if (eltype=='qu8' & rem(j,2) | rem(k,2))|sum(eltype~='qu8')  
                           % for eight node skipp midle
         % check that no lines vertical --------------------------------
         if X(ny,k)-X(1,k)~=0 
            k21=(Y(ny,k)-Y(1,k))/(X(ny,k)-X(1,k));
            m21=(Y(ny,k)*X(1,k)-X(ny,k)*Y(1,k))/(X(1,k)-X(ny,k));
         end
         if X(j,nx)-X(j,1)~=0
            k11=(Y(j,nx)-Y(j,1))/(X(j,nx)-X(j,1));
            m11=(Y(j,nx)*X(j,1)-X(j,nx)*Y(j,1))/(X(j,1)-X(j,nx));
         end
	 if X(ny,k)-X(1,k)==0  % line four vertical
            X(j,k)=X(1,k);
	    Y(j,k)=k11*X(1,k)+m11;
         elseif X(j,nx)-X(j,1)==0  % line two vertical
            X(j,k)=X(j,1);
	    Y(j,k)=k21*X(j,1)+m21;
         else
            X(j,k)=-(m11-m21)/(k11-k21);
	    Y(j,k)=(m11*k21-k11*m21)/(k21-k11);
         end
         nodes(j,k)=nextnode;
         CALFEM_COORD=[CALFEM_COORD; X(j,k) Y(j,k) 0];
         CALFEM_DOF=[CALFEM_DOF; (nextnode-1)*dofn+1:nextnode*dofn...
	        zeros(1,dfsb-dofn) ];
         nextnode=nextnode+1;
       end
      end
   end
end
% last line (number 3)
[x3,y3]=meshline(borders(3),dofn);
nodes(ny,2:nx-1)    =CALFEM_LINES_NODES(borders(3),3:nx);


