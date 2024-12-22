function err=linecombines(snr,lcomb,borders,eltype,dofn,npel)
% Treat linecombines in surface snr

global CALFEM_LINES
global CALFEM_LINES_NODES
global CALFEM_POINTS
global CALFEM_POINTS_NODE
global CALFEM_SURFS
global CALFEM_SURFS_ELEMS

if lcomb<3 
   if CALFEM_LINES(borders(lcomb+2),2)==0  % opposite line combined
      disp(['Sorry, geometry not able to treat opposite linecombines yet'])
      break;
   end
else   % to get linecombine right it has to be number 1 or 2
   borders=[borders(3) borders(4) borders(1) borders(2)];
   CALFEM_SURFS(snr,:)=[snr borders];
   lcomb=lcomb-2;
end

nextline=size(CALFEM_LINES,1)+1;
nextpoint=size(CALFEM_POINTS,1)+1;
nextsurf=size(CALFEM_SURFS,1)+1;

if CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),3)==CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),3)...
   | CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),4)==CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),4)
   CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),3:4)=...
           fliplr(CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),3:4));
   n=CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),2)+2;
   if size(CALFEM_LINES_NODES,2)>n-1     % lines have been meshed before
      CALFEM_LINES_NODES(CALFEM_LINES(borders(lcomb),4),2:n)=...
           fliplr(CALFEM_LINES_NODES(CALFEM_LINES(borders(lcomb),4),2:n));
   end
   CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),6)=-CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),6);
end

% make dummyline to direct boarders
if CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),3)~=CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),4)
   glines(nextline,0,CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),3),...
	CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),4),0,0,0);
else
   glines(nextline,0,CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),4),...
        CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),3),0,0,0);
end
CALFEM_LINES_NODES(nextline,1)=nextline;
newborders=borders;
newborders(lcomb)=nextline;

[flipp,err]=redirectlines(newborders);
if flipp(lcomb)
   CALFEM_LINES(borders(lcomb),3:4)=fliplr(CALFEM_LINES(borders(lcomb),3:4));
end
CALFEM_LINES(nextline,:)=[] ;  % take dummyline away
CALFEM_LINES_NODES(nextline,:)=[];   % take dummyline away

% check that combined line in right direction
if CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),3)==CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),4)
   CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),3:4)=...
           fliplr(CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),3:4));
   CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),6)=-CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),6);
   CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),3:4)=...
           fliplr(CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),3:4));
   CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),6)=-CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),6);
   if find(CALFEM_LINES_NODES(:,1)==CALFEM_LINES(borders(lcomb),3))   % if fliped line been meshed
      maxf=CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),2)+2;
      CALFEM_LINES_NODES(CALFEM_LINES(borders(lcomb),3),2:maxf)=...
   	   fliplr(CALFEM_LINES_NODES(CALFEM_LINES(borders(lcomb),3),2:maxf));
   end
   if find(CALFEM_LINES_NODES(:,1)==CALFEM_LINES(borders(lcomb),4))
      maxs=CALFEM_LINES(CALFEM_LINES(borders(lcomb),4),2)+2;
      CALFEM_LINES_NODES(CALFEM_LINES(borders(lcomb),4),2:maxs)=...
   	   fliplr(CALFEM_LINES_NODES(CALFEM_LINES(borders(lcomb),4),2:maxs));
   end
end

% create new point line and surf
opos=rem(lcomb+2,4);   % oposite line
lline=opos+1;
rline=rem(lcomb+1,4);
if ~opos opos=4; end
if ~rline rline=4; end
divy=CALFEM_LINES(borders(lline),2);

divfirst=CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),2);
[xop,yop]=divline(borders(opos));
gpoints(nextpoint,xop(divfirst+1),yop(divfirst+1),0,0);
% two new lines on opposit side
glines(nextline,divfirst,CALFEM_LINES(borders(opos),3),nextpoint,CALFEM_LINES(borders(opos),5),0,0);
glines(nextline+1,CALFEM_LINES(borders(opos),2)-divfirst,nextpoint,CALFEM_LINES(borders(opos),4),...
	CALFEM_LINES(borders(opos),5),0,0);
% line dividing into two surfaces
glines(nextline+2,divy,CALFEM_LINES(CALFEM_LINES(borders(lcomb),3),4),nextpoint,0,0,0);
% two new surfaces
gsurfs(nextsurf,CALFEM_LINES(borders(lcomb),3),nextline+2,nextline,borders(lline));
gsurfs(nextsurf+1,CALFEM_LINES(borders(lcomb),4),borders(rline),nextline+1,nextline+2);

% mesh the two new surfaces
err=meshsurf(nextsurf,eltype,dofn,npel);
err=meshsurf(nextsurf+1,eltype,dofn,npel);
if err break; end

% clean up
[nlns,mnpl]=size(CALFEM_LINES_NODES);
newline=[borders(opos) CALFEM_LINES_NODES(nextline,2:divfirst+1) ... 
             CALFEM_LINES_NODES(nextline+1,2:CALFEM_LINES(borders(opos),2)-divfirst+2)];
nnpl=size(newline,2);
if nnpl>mnpl
   CALFEM_LINES_NODES=[CALFEM_LINES_NODES zeros(nlns,nnpl-mnpl)];
else
   newline=[newline zeros(1,mnpl-nnpl)];
end
CALFEM_LINES_NODES(borders(opos),:)=newline;
CALFEM_POINTS(nextpoint,:)=[];
CALFEM_POINTS_NODE(nextpoint,:)=[];
CALFEM_LINES(nextline:nextline+2,:)=[];
CALFEM_LINES_NODES(nextline:nextline+2,:)=[];
CALFEM_SURFS(nextsurf:nextsurf+1,:)=[];
lnr=find(CALFEM_SURFS_ELEMS(:,1)==nextsurf);
frst=min([CALFEM_SURFS_ELEMS(lnr,2:3) CALFEM_SURFS_ELEMS(lnr+1,2:3)]);
lst=max([CALFEM_SURFS_ELEMS(lnr,2:3) CALFEM_SURFS_ELEMS(lnr+1,2:3)]);
CALFEM_SURFS_ELEMS=[CALFEM_SURFS_ELEMS; snr frst lst];
CALFEM_SURFS_ELEMS(lnr:lnr+1,:)=[];
