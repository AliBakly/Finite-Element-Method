function [s1,s2]=pltstyle(plotpar)
%-------------------------------------------------------------
% PURPOSE 
%   Define define linetype,linecolor and markertype character codes. 
%
% INPUT 
%    plotpar=[ linetype, linecolor, nodemark ]
% 
%             linetype=1 -> solid    linecolor=1 -> white
%                      2 -> dashed             2 -> green
%                      3 -> dotted             3 -> yellow
%                                              4 -> red
%                                              5 -> blue
%                                              6 -> black
%                                              7 -> magenta
%
%             nodemark=1 -> circle       
%                      2 -> star               
%                      0 -> no mark             
%                      3 -> no mark at all             
% OUTPUT
%     s1: linetype and color for mesh lines
%     s2: type and color for node markers
%-------------------------------------------------------------

% LAST MODIFIED: P-E Austrell 1993-12-14 
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if plotpar(1)==1 ; s1='-';
 elseif plotpar(1)==2 ; s1='--';
 elseif plotpar(1)==3 ; s1=':';
 else error('??? Error in variable plotpar(1)!');
      %break;
 end
 
 if plotpar(2)==1 ; s1=[s1,'w'];
 elseif plotpar(2)==2 ; s1=[s1,'g'];
 elseif plotpar(2)==3 ; s1=[s1,'y'];
 elseif plotpar(2)==4 ; s1=[s1,'r'];
 elseif plotpar(2)==5 ; s1=[s1,'b'];
 elseif plotpar(2)==6 ; s1=[s1,'k'];
 elseif plotpar(2)==7 ; s1=[s1,'m'];
else error('??? Error in variable plotpar(2)!');
      
 end
 
 if plotpar(3)==1 ; s2='wo';
 elseif plotpar(3)==2 ; s2='w*';
 elseif plotpar(3)==0 ; s2='w.';
 elseif plotpar(3)==3 ; s2=' ';
 else error('??? Error in variable plotpar(3)!');
      
 end
%--------------------------end-------------------------------- 
 
