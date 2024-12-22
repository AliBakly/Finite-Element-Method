function [flipp,err]=redirectlines(borders)
% make shure lines in right directions, borderline 1
% from line 4 to line 2, line 3 in the same direction 
% and line 2 and 4 from line 1 to 3.
% if err = 1 the lines do not form a closed surface
% if err = 2 the lines orderd clockwise (wrong way)
global CALFEM_POINTS
global CALFEM_LINES
global CALFEM_LINES_NODES

if isempty(CALFEM_LINES_NODES)
   CALFEM_LINES_NODES=0;
end
p1=CALFEM_LINES(borders(1),3);
p2=CALFEM_LINES(borders(1),4);
p3=CALFEM_LINES(borders(2),3);
p4=CALFEM_LINES(borders(2),4);
p5=CALFEM_LINES(borders(3),3);
p6=CALFEM_LINES(borders(3),4);
p7=CALFEM_LINES(borders(4),3);
p8=CALFEM_LINES(borders(4),4);

err=0;
flipp=zeros(4,1);

if p2==p7 | p2==p8     % line 1 wrong way
   CALFEM_LINES(borders(1),3:4)=fliplr(CALFEM_LINES(borders(1),3:4));
   CALFEM_LINES(borders(1),6)=-CALFEM_LINES(borders(1),6);
   p1=p2; p2=CALFEM_LINES(borders(1),4); 
   if find(CALFEM_LINES_NODES(:,1)==borders(1))
      flipp(1)=1;  % for linecombines
      n=CALFEM_LINES(borders(1),2)+2;
      if size(CALFEM_LINES_NODES,2)>n-1     % realy meshed before
         CALFEM_LINES_NODES(borders(1),2:n)=fliplr(CALFEM_LINES_NODES(borders(1),2:n));
      end
   end
end
if p2==p4        % line 2 wrong way
   CALFEM_LINES(borders(2),3:4)=fliplr(CALFEM_LINES(borders(2),3:4));
   CALFEM_LINES(borders(2),6)=-CALFEM_LINES(borders(2),6);
   p3=p4; p4=CALFEM_LINES(borders(2),4);
   if find(CALFEM_LINES_NODES(:,1)==borders(2))
      flipp(2)=1;
      n=CALFEM_LINES(borders(2),2)+2;
      if size(CALFEM_LINES_NODES,2)>=n-1
	 CALFEM_LINES_NODES(borders(2),2:n)=fliplr(CALFEM_LINES_NODES(borders(2),2:n));
      end
   end
end
if p4==p5	 % line 3 wrong way
   CALFEM_LINES(borders(3),3:4)=fliplr(CALFEM_LINES(borders(3),3:4));
   CALFEM_LINES(borders(3),6)=-CALFEM_LINES(borders(3),6);
   p5=p6; p6=CALFEM_LINES(borders(3),4);
   if find(CALFEM_LINES_NODES(:,1)==borders(3))
      flipp(3)=1;  
      n=CALFEM_LINES(borders(3),2)+2;
      if size(CALFEM_LINES_NODES,2)>n-1
	 CALFEM_LINES_NODES(borders(3),2:n)=fliplr(CALFEM_LINES_NODES(borders(3),2:n));
      end
   end
end
if p8==p1	 % line 4 wrong way
   CALFEM_LINES(borders(4),3:4)=fliplr(CALFEM_LINES(borders(4),3:4));
   CALFEM_LINES(borders(4),6)=-CALFEM_LINES(borders(4),6);
   p7=p8; p8=CALFEM_LINES(borders(4),4);
   if find(CALFEM_LINES_NODES(:,1)==borders(4))
      flipp(4)=1;  
      n=CALFEM_LINES(borders(4),2)+2;
      if size(CALFEM_LINES_NODES,2)>n-1
	 CALFEM_LINES_NODES(borders(4),2:n)=fliplr(CALFEM_LINES_NODES(borders(4),2:n));
      end
   end
end

if p2~=p3 | p4~=p6 | p5~=p8 | p7~=p1    % not closed
   err=1; break
end
ec=[CALFEM_POINTS(p1,2) CALFEM_POINTS(p3,2) CALFEM_POINTS(p6,2) CALFEM_POINTS(p8,2);
    CALFEM_POINTS(p1,3) CALFEM_POINTS(p3,3) CALFEM_POINTS(p6,3) CALFEM_POINTS(p8,3)];
[gp,wp]=gauss(1,2);
[N,dNr,blk]=shape(4,gp);
JT=dNr*ec';
detJ=det(JT(blk(:,1),:));
if detJ<0
   err=2;
end