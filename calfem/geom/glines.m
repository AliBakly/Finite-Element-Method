function l=glines(lnum,div,p1,p2,p3,rel,dr)
% l=glines(lnum,div,p1,p2,p3,rel,dr)
%------------------------------------------------------------------------------
% PURPOSE
%  Create a line which is used to define the geometry
%
% INPUT:  pnum     number of line, if 0 a number will be generated
%
%         div      number of elements to divide the line
%
%         p1,p2,p3 points to define start, end [and center] of line
%
%         rel      relation between largest and smalles element size
%
%         dr       if the line is not to be drawn dr=0
%
% OUTPUT: l        if no arguments the list of lines will be returned
%-------------------------------------------------------------

% REFERENCES
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology

global CALFEM_LINES;
if nargin==0
   l=CALFEM_LINES;
   return
end
if div==0 div=8; end
if nargin==4 p3=0; rel=0; dr=1; end
if nargin==5 rel=0; dr=1; end
if nargin==6 dr=1; end
nrlines=size(p1,1);
if nrlines>1     % many lines specified in a matrix
   for i=1:nrlines
      if size(lnum,1)<i lnum(i)=0; end
      if size(div,1)<i div(i)=0; end
      if size(p3,1)<i p3(i)=0; end
      if size(rel,1)<i rel(i)=0; end
      if size(dr,1)<i dr(i)=1; end
      glines(lnum(i),div(i),p1(i),p2(i),p3(i),rel(i),dr(i));
   end
   break
end
[r,c]=size(CALFEM_LINES);
if r==0
   if lnum==0
      lnum=1;
   end
   CALFEM_LINES(lnum,1:6)=[lnum div p1 p2 p3 rel]; 
   if dr gdraw('lin',1,lnum); end
else
   if lnum==0
      n=max(CALFEM_LINES(:,1));
      lnum=n+1;
      disp(['Linenumber ',num2str(lnum),' created'])
   end
   if isempty(find(CALFEM_LINES(:,1)==lnum))
      CALFEM_LINES(lnum,:)=[lnum div p1 p2 p3 rel];
      if dr gdraw('lin',1,lnum); end
   else
      disp('GEOMETRY ERROR: linenumber already in use')
   end
end
