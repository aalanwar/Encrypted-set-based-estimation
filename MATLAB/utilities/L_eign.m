function [x] = L_eign(fstate,x,hmeas,y)


[x1,F]=jaccsd(fstate,x);    %nonlinear update and linearization at current state
[z1,H]=jaccsd(hmeas,x1);    %nonlinear measurement and linearization
% L0=ones(3,1);
% options = optimoptions(@fminunc,'Algorithm', 'quasi-newton','Display','off');
% %find the weights
% L = fminunc(@fun,L0, options);


L0=0.1*ones(3,1);
A = [];
b = [];
Aeq = [];
beq = [];
%lb = -6*ones(1,length(x0));
%ub = 6*ones(1,length(x0));
lb = [];
ub = [];
nonlcon = @circlecon;
options = optimoptions(@fmincon,'Algorithm', 'sqp', 'TolX', 1e-9, 'TolFun', 1e-9,'Display','off');
L = fmincon(@fun,L0,A,b,Aeq,beq,lb,ub,nonlcon, options);


x = F*x + L *(y-H*x);


    function nfro = fun(L)
        
        %syms z;
        %numerator = H;
        %denominator = eye(3)*z+-F+L*H;
        %ts = 0.1;
        %sys = tf(numerator,denominator,ts);
        ts = 0.01;
        if H(1) ==1
            %H 1 0 0
            tf_n1 = tf(1,[1 L(1)-1],ts);
            tf_n2 = tf(0,1,ts);
            tf_n3 = tf(0,1,ts);
            
            tf_v = tf(L(1),[1 L(1)-1],ts);
        elseif H(2)==1
            tf_n1 = tf(0,1,ts);
            tf_n2 = tf(1,[1 L(2)-1],ts);
            tf_n3 = tf(0,1,ts);
            
            tf_v = tf(L(2),[1 L(2)-1],ts);
        else
            tf_n1 = tf(0,1,ts);
            tf_n2 = tf(0,1,ts);
            tf_n3 = tf(1,[1 L(3)-1],ts)  ;
            
            tf_v = tf(L(3),[1 L(3)-1],ts);
        end
        tf_p1 = [tf_n1,tf_n2,tf_n3];
        nfro_p1 = norm(tf_p1,Inf);
        
        tf_p2 = 1 - tf_v;
        nfro_p2 = norm(tf_p2,Inf);
        %1 - L2/(L2 + z - 1)
        nfro = nfro_p1*nfro_p2;
    end
    function [c,ceq] = circlecon(L)
        % c is less than zero c<=0
        eig_value= eig(F-L*H);
        c=norm(eig_value - ones(length(eig_value),1));
        
        
        ceq =[];%sum(lambda)-1;
        
    end
end

function [z,A]=jaccsd(fun,x)
% JACCSD Jacobian through complex step differentiation
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
z=fun(x);
n=numel(x);
m=numel(z);
A=zeros(m,n);
h=n*eps;
for k=1:n
    x1=x;
    x1(k)=x1(k)+h*1i;
    A(:,k)=imag(fun(x1))/h;
end
end