function err=meshsurf(snr,eltype,dofn,npel)
% assign nodenumbers for a surface

global CALFEM_DOF
global CALFEM_COORD
global CALFEM_EDOF
global CALFEM_LINES
global CALFEM_LINES_NODES
global CALFEM_SURFS
global CALFEM_SURFS_ELEMS

borders=CALFEM_SURFS(snr,2:5);  % border lines
err=0;
% check for line-combines
lc=0;
for k=1:4
   if CALFEM_LINES(borders(k),2)==0 
      err=linecombines(snr,k,borders,eltype,dofn,npel); lc=1;break;
   end
end

if lc | err break; end

[flipp,err]=redirectlines(borders);	% make sure lines in right direction

if err==1
   disp(['Lines do not form a closed surface for surface nr ',num2str(snr)])
   break
elseif err==2
   % borders was in the wrong way
   CALFEM_SURFS(snr,2:5)=fliplr(CALFEM_SURFS(snr,2:5));
   err=meshsurf(snr,eltype,dofn,npel);
   break
end

[X,Y,nodes]=surfcoords(borders,eltype,dofn);  % coordinates for nodes in surface

nx=CALFEM_LINES(borders(1),2); ny=CALFEM_LINES(borders(2),2);
[nextelement,maxdofs]=size(CALFEM_EDOF);
nextelement=nextelement+1;
firstelement=nextelement;
maxdofs=maxdofs-1;
thisdofs=dofn*npel;
if thisdofs>maxdofs      % more dofs than before
   CALFEM_EDOF=[CALFEM_EDOF zeros(nextelement-1,thisdofs-maxdofs)];
end
sros=zeros(1,maxdofs-thisdofs);
step=1;
if eltype=='qu8'|eltype=='qr6' 
   step=2;   % each element takes two divides
end
nextdof=max(max(CALFEM_DOF))+1;
Dof=CALFEM_DOF;
for r=1:step:ny
   for k=1:step:nx
      if eltype=='qu4'
	 CALFEM_EDOF(nextelement,:)=[nextelement ...
	      ((nodes(r,k)-1)*dofn+1):nodes(r,k)*dofn ...
	      ((nodes(r,k+1)-1)*dofn+1):nodes(r,k+1)*dofn ...
	      ((nodes(r+1,k+1)-1)*dofn+1):nodes(r+1,k+1)*dofn ...
	      ((nodes(r+1,k)-1)*dofn+1):nodes(r+1,k)*dofn sros];
         nextelement=nextelement+1;
      elseif eltype=='tr3'
         CALFEM_EDOF(nextelement,:)=[nextelement ...
              ((nodes(r,k)-1)*dofn+1):nodes(r,k)*dofn ...
              ((nodes(r,k+1)-1)*dofn+1):nodes(r,k+1)*dofn ...
              ((nodes(r+1,k+1)-1)*dofn+1):nodes(r+1,k+1)*dofn sros];
         CALFEM_EDOF(nextelement+1,:)=[nextelement+1 ...
              ((nodes(r,k)-1)*dofn+1):nodes(r,k)*dofn ...
              ((nodes(r+1,k+1)-1)*dofn+1):nodes(r+1,k+1)*dofn ...
              ((nodes(r+1,k)-1)*dofn+1):nodes(r+1,k)*dofn sros];
	 nextelement=nextelement+2;
      elseif eltype=='qu8'
         CALFEM_EDOF(nextelement,:)=[nextelement ...
	      ((nodes(r,k)-1)*dofn+1):nodes(r,k)*dofn...
	      ((nodes(r,k+2)-1)*dofn+1):nodes(r,k+2)*dofn...
	      ((nodes(r+2,k+2)-1)*dofn+1):nodes(r+2,k+2)*dofn...
	      ((nodes(r+2,k)-1)*dofn+1):nodes(r+2,k)*dofn...
	      ((nodes(r,k+1)-1)*dofn+1):nodes(r,k+1)*dofn...
	      ((nodes(r+1,k+2)-1)*dofn+1):nodes(r+1,k+2)*dofn...
	      ((nodes(r+2,k+1)-1)*dofn+1):nodes(r+2,k+1)*dofn...
	      ((nodes(r+1,k)-1)*dofn+1):nodes(r+1,k)*dofn sros];
	 nextelement=nextelement+1;
      elseif eltype=='qr3'
         cent=(CALFEM_COORD(nodes(r,k),:)+CALFEM_COORD(nodes(r,k+1),:)+...
               CALFEM_COORD(nodes(r+1,k+1),:)+CALFEM_COORD(nodes(r+1,k),:))/4;
	 CALFEM_COORD=[CALFEM_COORD;cent];
	 nnode=nextdof:nextdof+dofn-1;
	 CALFEM_DOF=[CALFEM_DOF; nnode];
	 CALFEM_EDOF(nextelement,:)=[nextelement ...
	      Dof(nodes(r,k),:) Dof(nodes(r,k+1),:) nnode sros];
	 CALFEM_EDOF(nextelement+1,:)=[nextelement+1 ...
	      Dof(nodes(r,k+1),:) Dof(nodes(r+1,k+1),:) nnode sros];
	 CALFEM_EDOF(nextelement+2,:)=[nextelement+2 ...
	      Dof(nodes(r+1,k+1),:) Dof(nodes(r+1,k),:) nnode sros];
	 CALFEM_EDOF(nextelement+3,:)=[nextelement+3 ...
	      Dof(nodes(r+1,k),:) Dof(nodes(r,k),:) nnode sros];
	 nextelement=nextelement+4;
	 nextdof=nextdof+dofn;
      elseif eltype=='qr6'
         c1=(CALFEM_COORD(nodes(r,k),:)+CALFEM_COORD(nodes(r+1,k+1),:))/2;
         c2=(CALFEM_COORD(nodes(r,k+2),:)+CALFEM_COORD(nodes(r+1,k+1),:))/2;
         c3=(CALFEM_COORD(nodes(r+2,k+2),:)+CALFEM_COORD(nodes(r+1,k+1),:))/2;
         c4=(CALFEM_COORD(nodes(r+2,k),:)+CALFEM_COORD(nodes(r+1,k+1),:))/2;
         CALFEM_COORD=[CALFEM_COORD;c1;c2;c3;c4];
	 nnode1=nextdof:nextdof+dofn-1;
	 nnode2=nnode1+dofn; nnode3=nnode2+dofn; nnode4=nnode3+dofn;
	 CALFEM_DOF=[CALFEM_DOF;nnode1;nnode2;nnode3;nnode4];
	 CALFEM_EDOF(nextelement,:)=[nextelement ...
              Dof(nodes(r,k),:)      Dof(nodes(r,k+2),:)   Dof(nodes(r+1,k+1),:)...
	      Dof(nodes(r,k+1),:)    nnode2                nnode1        sros];
	 CALFEM_EDOF(nextelement+1,:)=[nextelement+1 ...
              Dof(nodes(r,k+2),:)    Dof(nodes(r+2,k+2),:) Dof(nodes(r+1,k+1),:)...
	      Dof(nodes(r+1,k+2),:)  nnode3                nnode2        sros];
	 CALFEM_EDOF(nextelement+2,:)=[nextelement+2 ...
              Dof(nodes(r+2,k+2),:)  Dof(nodes(r+2,k),:)   Dof(nodes(r+1,k+1),:)...
	      Dof(nodes(r+2,k+1),:)  nnode4                nnode3        sros];
	 CALFEM_EDOF(nextelement+3,:)=[nextelement+3 ...
              Dof(nodes(r+2,k),:)    Dof(nodes(r,k),:)     Dof(nodes(r+1,k+1),:)...
	      Dof(nodes(r+1,k),:)    nnode1                nnode4        sros];
	 nextelement=nextelement+4;
	 nextdof=nextdof+dofn*4;
      end
   end
end
CALFEM_SURFS_ELEMS=[CALFEM_SURFS_ELEMS; snr firstelement nextelement-1];