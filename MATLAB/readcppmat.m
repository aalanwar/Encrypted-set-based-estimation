clear all
load('CMatFiles/ZonoStrips.mat')
%get the number of var in workspace
variablesInCurrentWorkspace = who;
numVariables = length(variablesInCurrentWorkspace);
index =1;
numoffiles = length(dir(fullfile('MatFiles/', '*.mat')));
for i=0:numVariables/2-1
 c_i= eval(strcat('c_',num2str(i)));
 G_i= eval(strcat('G_',num2str(i)));
 zonol{index} = zonotope([c_i,G_i]);
 sup(:,index) = supremum(interval( zonol{index}));
 infi(:,index) = infimum(interval( zonol{index}));
 cen (:,index) = c_i;
 index = index +1;
end



figure; 
plot(cen(1,:))
hold on 
plot(cen(2,:))
plot(cen(3,:))