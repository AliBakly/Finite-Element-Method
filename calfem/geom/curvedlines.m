function [X,Y,nodes]=curvedlines(cul,borders,X,Y,nodes,dofn)
% make mesh for curved lines (internal nodes)

if cul(1) & cul(3)    % line 1 and 3 curved
   close=1; oppos=3; left=4; right=2;
   [X,Y,nodes]=twocurved(X,Y,nodes,borders,dofn,close,oppos,left,right);
elseif cul(2) & cul(4)    % line 2 and 4 curved
   close=2; oppos=4; left=1; right=3;
   [X,Y,nodes]=twocurved(X',Y',nodes',borders,dofn,close,oppos,left,right);
   X=X'; Y=Y'; nodes=nodes';
elseif cul(1)   
   close=1; oppos=3; left=4; right=2;   
   [X,Y,nodes]=onecurved(X,Y,nodes,borders,dofn,close,oppos,left,right);
elseif cul(2)   
   close=2; oppos=4; left=1; right=3;   
   [X,Y,nodes]=onecurved(X',Y',nodes',borders,dofn,close,oppos,left,right);
   X=X'; Y=Y'; nodes=nodes';
elseif cul(3)   
   close=3; oppos=1; left=2; right=4;   
   [X,Y,nodes]=onecurved(X,Y,nodes,borders,dofn,close,oppos,left,right);
elseif cul(4)   
   close=4; oppos=2; left=3; right=1;   
   [X,Y,nodes]=onecurved(X',Y',nodes',borders,dofn,close,oppos,left,right);
   X=X'; Y=Y'; nodes=nodes';
end
