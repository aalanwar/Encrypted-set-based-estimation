clear all
load('<path-to-mat-file');
n = size(cc,1);
m = size(GG,2);
N = size(cc,2);
close all
figure
hold on
for i=1:N
    G = GG(1+(i-1)*n:i*n,:);
    c = cc(:,i);
    x = xx(:,i);
    plot(zonotope([c,G]), [1,2], 'k');
    plot(x(1),x(2),'xk');
    plot(0,0,'or');
end