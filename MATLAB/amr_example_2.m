clear all
load('<path-to-mat-file');
n = size(cc,1);
m = size(GG,2);
N = size(cc,2)/(e+1);
close all
figure
hold on
plot(0,0,'or');
for i=1:N
    G = GG(1+(i-1)*(e+1)*n:(i-1)*(e+1)*n + n,:);
    c = cc(:,1+(i-1)*(e+1));
    x = xx(:,i);
    plot(x(1),x(2),'xk');
    H = [];
    K = [];
    ZZ = {};
    for j=1:e
        Gint = GG(1+(i-1)*(e+1)*n + j*n:(i-1)*(e+1)*n + j*n + n,:);
        cint = cc(:,1+(i-1)*(e+1)+j);
        Ztmp = zonotope([cint,Gint]);
        ZZ{end+1} = Ztmp;
        HS = halfspace(Ztmp);
        H = [H;HS.halfspace.H];
        K = [K;HS.halfspace.K];
        plot(Ztmp);
    end
    P = Polyhedron(H,K);
    P.plot('wire',1,'LineStyle','--');
    plot(zonotope([c,G]), [1,2], 'Color','r','LineWidth',2);
    %%compute intersection in MATLAB
    A = zeros(e,e);
    for j=1:e
        Gj = ZZ{j}.Z(:,2:end);
        for k=j:e
            Gk = ZZ{k}.Z(:,2:end);
            tmp = trace(Gj*Gk');
            A(j,k) = tmp;
            A(k,j) = tmp;
        end
    end
    w = inv(A)*ones(e,1)/(ones(1,e)*inv(A)*ones(e,1));
    G_m = [];
    c_m = zeros(n,1);
    for j=1:length(ZZ)
        G_m = [G_m,w(j)/sum(w)*ZZ{j}.Z(:,2:end)];
        c_m = c_m + w(j)/sum(w)*ZZ{j}.Z(:,1);
    end
    plot(zonotope([c_m,G_m]),[1,2],'Color','g','LineWidth',2,'LineStyle','--');
end