function [res_conzonotope]=intersectConZonoStrip1(z1,hl,Rl,yl)
% intersectZonoStrip - computes the intersection between one constrained zonotope and list of strips 
%
%% the strip is defined as | hx-y | <= d
%% example with three strips and one zonotope:
% hl{1} = [1 0];
% Rl{1} = 5;
% yl{1} = -2;
% 
% hl{2} = [0 1];
% Rl{2} = 3;
% yl{2} = 2;
% 
% hl{3} = [1 1];
% Rl{3} = 3;
% yl{3} = 2;
% 
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%   cZono1 = conZonotope(Z,A,b);
%res_czono= untitled(cZono1,hl,Rl,yl);

%%just for plotting strips
%poly = mptPolytope([1 0;-1 0; 0 1;0 -1; 1 1;-1 -1],[3;7;5;1;5;1]);


%figure; hold on 
%plot(cZono1,[1 2],'r-+');
%plot(poly,[1 2],'r-*');
%plot(res_czono,[1 2],'b-*');

%legend('zonotope','strips','czonoStrips');
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author: Amr Alanwar
% Written: 9-Mar-2020
% Last update: ---
%              
% Last revision: ---

%------------- BEGIN CODE --------------


lambda = 0.5*zeros(length(z1.center),length(Rl));
%prepare center
c_new=z1.Z(:,1);
for i=1:length( Rl)
    c_new = c_new + lambda(:,i)*( yl{i} - hl{i}*z1.Z(:,1) );
end

%prepare generators
part1 = eye(length(z1.Z(:,1)));
if isempty(z1.A)
    A_new =[];
    b_new =[];
else
    A_new = [ z1.A zeros(size(z1.A,1),length( Rl))];
    b_new = z1.b;
end

for ii=1:length(Rl)
    part1 = part1 - lambda(:,ii)*hl{ii};
    part2(:,ii) = Rl{ii}*lambda(:,ii);
    A_new = [A_new ; hl{ii}*z1.Z(:,2:end) , zeros(1,ii-1),-Rl{ii},zeros(1,length(Rl)-ii)];
    b_new = [b_new; yl{ii}-(hl{ii}*z1.Z(:,1))];
end
part1 = part1 * z1.Z(:,2:end);
H_new = [part1 part2];

res_conzonotope = conZonotope([c_new H_new],A_new,b_new);





end
