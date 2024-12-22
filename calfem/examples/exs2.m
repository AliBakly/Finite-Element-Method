% example exs2 
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of one dimensional heat flow.
%----------------------------------------------------------------

% REFERENCES
%     P-E Austrell 1994-03-08 
%     K-G Olsson 1995-09-28
%----------------------------------------------------------------
echo on
 
%----- Topology matrix Edof -------------------------------------

 Edof=[1 1 2;
       2 2 3;
       3 3 4;
       4 4 5;
       5 5 6];
 
%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(6); 
 f=zeros(6,1);

%----- Element properties ---------------------------------------

 ep1=[ 1/0.07 ];	ep2=[ 1.7/0.07 ];
 ep3=[ 0.040/0.10 ];	ep4=[ 1.7/0.10 ];
 ep5=[ 1/0.18 ];
 
%----- Element stiffness matrices  ------------------------------

 Ke1=spring1e(ep1);	Ke2=spring1e(ep2);
 Ke3=spring1e(ep3);	Ke4=spring1e(ep4);
 Ke5=spring1e(ep5);

 
%----- Assemble Ke into K ---------------------------------------

 K=assem(Edof(1,:),K,Ke1);	K=assem(Edof(2,:),K,Ke2); 
 K=assem(Edof(3,:),K,Ke3);	K=assem(Edof(4,:),K,Ke4);
 K=assem(Edof(5,:),K,Ke5);
 
%----- Solve the system of equations ----------------------------

 bc= [1 -17; 6 20];   
 T=solveq(K,f,bc)

%----- Element flows -------------------------------------------

 ed1=extract(Edof(1,:),T);
 ed2=extract(Edof(2,:),T);
 ed3=extract(Edof(3,:),T);
 ed4=extract(Edof(4,:),T);
 ed5=extract(Edof(5,:),T);

 q1=spring1s(ep1,ed1)
 q2=spring1s(ep2,ed2)
 q3=spring1s(ep3,ed3)
 q4=spring1s(ep4,ed4)
 q5=spring1s(ep5,ed5)

%---------------------------- end -------------------------------
 echo off
