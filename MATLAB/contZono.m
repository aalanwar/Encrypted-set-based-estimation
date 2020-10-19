   
clear all
% constrained zonotopes
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%    cZono1 = conZonotope(Z,A,b);

%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1];
%    b = 1;
%    cZono2 = conZonotope(Z,A,b);
% 
% res1 = Z & cZono2;

% hl{1} = [1 0];
% Rl{1} = 5;
% yl{1} = -2;
% 
% hl{2} = [0 1];
% Rl{2} = 3;
% yl{2} = 2;
% 
% res_conzono= intersectConZonoStrip1(cZono1,hl,Rl,yl);

%  redZono = reduce(cZono1,'girard',1,0);
 
Z = [0 3 0 1;0 0 2 1; 1 0 0 0];
A = [1 0 1];
b = 1;
cZono1 = conZonotope(Z,A,b);

newceter = [ Z(:,1); -b];
newgen = [ Z(:,2:end); A];
liftzono= zonotope([newceter,newgen]);

figure;plot(cZono1)
hold on ;plot(liftzono,[1 2],'r')

Q2 = diag([1 1 1]);
cQ2 = [ 0;0;0];
A = [0 0 0];
b = 0;
Qconzono = conZonotope([cQ2 Q2],A,b);

timeup = cZono1 + Qconzono;
figure;plot(cZono1)
hold on ;plot(timeup,[1 2],'r')


center  =zeros(3,1);
width   =7.5;%16
x_conzono=conZonotope(zonotope([center,width*eye(3,3)]));
