function l=glinec(lnum,l1,l2,l3)
% l=glinec(lnum,l1,l2,l3)
%------------------------------------------------------------------------------
% PURPOSE
%  Define a line from already created lines
%
% INPUT:  lnum     number of line, if 0 a number will be generated
%
%         l1,l2,l3 lines to combine
%
% OUTPUT: l        if no arguments the list of lines will be returned
%-------------------------------------------------------------

% REFERENCES
%   M Ristinmaa 1994-05-24
%   A Magnusson 1994-05-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
global CALFEM_LINES;
if nargin==0
   l=CALFEM_LINES;
   return
end
if nargin==3 l3=-1; end
nrlines=size(l1,1);
if nrlines>1     % many lines specified in a matrix
   for i=1:nrlines
      if size(lnum,1)<i lnum(i)=0; end
      if size(l3,1)<i l3(i)=0; end
      glinec(lnum(i),l1(i),l2(i),l3(i));
   end
   break
end
test=[find(CALFEM_LINES(:,1)==l1) find(CALFEM_LINES(:,1)==l2)
      find(CALFEM_LINES(:,1)==l3)];
if size(test)~=[1 3]
     if size(test)~=[1 2]
       disp('LINE ERROR: could not find line')
       return
     end
end 
if nargin==3 l3=0; end
[r,c]=size(CALFEM_LINES);
if lnum==0
  n=max(CALFEM_LINES(:,1));
  lnum=n+1;
  disp(['Linenumber ',num2str(lnum),' created'])
end
if isempty(find(CALFEM_LINES(:,1)==lnum))
  CALFEM_LINES(lnum,:)=[lnum 0 l1 l2 l3 0];
%  gdraw('lin',1,lnum);
else
  disp('GEOMETRY ERROR: linenumber already in use')
end
