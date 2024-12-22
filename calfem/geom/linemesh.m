function err=linemesh(lnr,eltype,dofn,npel)
% mesh a line

global CALFEM_EDOF
global CALFEM_LINES
global CALFEM_LINES_NODES
global CALFEM_LINES_ELEMS

err=0;
nextelem=size(CALFEM_EDOF,1)+1;
firstelem=nextelem;
div=CALFEM_LINES(lnr,2);
if div		% not a line-combine
   [x,y]=meshline(lnr,dofn);
   if npel==2
     for k=1:div
      CALFEM_EDOF(nextelem,:)=[nextelem ...
      (CALFEM_LINES_NODES(lnr,k+1)-1)*dofn+1:CALFEM_LINES_NODES(lnr,k+1)*dofn...
      (CALFEM_LINES_NODES(lnr,k+2)-1)*dofn+1:CALFEM_LINES_NODES(lnr,k+2)*dofn];
      nextelem=nextelem+1;
     end     
 
   elseif npel==3
     CALFEM_LINES_NODES
     for k=1:2:div
      CALFEM_EDOF(nextelem,:)=[nextelem ...
      (CALFEM_LINES_NODES(lnr,k+1)-1)*dofn+1:CALFEM_LINES_NODES(lnr,k+1)*dofn...
      (CALFEM_LINES_NODES(lnr,k+2)-1)*dofn+1:CALFEM_LINES_NODES(lnr,k+2)*dofn...
      (CALFEM_LINES_NODES(lnr,k+3)-1)*dofn+1:CALFEM_LINES_NODES(lnr,k+3)*dofn];
       nextelem=nextelem+1;
     end     
   end
else            % line-combine
   if CALFEM_LINES(CALFEM_LINES(lnr,3),3)==CALFEM_LINES(CALFEM_LINES(lnr,4),3)|...
      CALFEM_LINES(CALFEM_LINES(lnr,3),4)==CALFEM_LINES(CALFEM_LINES(lnr,4),4)
      % the lines not in same direction
      CALFEM_LINES(CALFEM_LINES(lnr,4),3:4)=...
                 fliplr(CALFEM_LINES(CALFEM_LINES(lnr,4),3:4));
      CALFEM_LINES(CALFEM_LINES(lnr,4),6)=-CALFEM_LINES(CALFEM_LINES(lnr,4),6);
   end
   gmesh('line',CALFEM_LINES(lnr,3),eltype,dofn);
   gmesh('line',CALFEM_LINES(lnr,4),eltype,dofn);
   nextelem=size(CALFEM_EDOF,1)+1;
end
CALFEM_LINES_ELEMS=[CALFEM_LINES_ELEMS; lnr firstelem nextelem-1];
