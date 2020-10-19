

%% ------ Question 1 ------- %%

clear all
close all
load('cache/conzonotest.mat')

figure ;plot(x_conzono,[1 2],'r')
xlabel('z1')
ylabel('z2')

figure ;plot(x_conzono,[1 3],'r')
xlabel('z1')
ylabel('z3')
x_conzono.center



%% ------ Question 2 ------- %%

 Q = diag([1 1 1]);
 cQ = [ 0;0;0];
 %   A = [0 0 0];
 %  b = 0;
 Qconzono = conZonotope(zonotope([cQ Q]));
 %Qconzono = conZonotope([cQ2 Q2],A,b);
 
 
 G = diag([2 2 2]);
 c = [ 2;1;3];
    A = [1 0 1];
   b = 1;
 conz1 = conZonotope([c G],A,b);
 
 res = conz1 + Qconzono