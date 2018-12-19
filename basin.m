close all
clear all
%
l = [-2,-1.5];
u = [2,1.5];
pf = [-.5,0;.5,0];
N = 50;%
T = 100;

[X,Y,B] = basin2D(l,u,pf,N,T);
%
[Xm,Ym] = meshgrid(X,Y);
%Bm = griddata(X,Y,B,Xm,Ym,'natural');

figure()
surf(Xm,Ym,B,'EdgeColor','None')
view(2)
%
function [X,Y,B] = basin2D(l,u,pf,N,T)
    f = waitbar(0,'Please wait...');
    X = linspace(l(1),u(1),N);
    Y = linspace(l(2),u(2),N);
    %
    Npf = size(pf,1);
    %
    B = zeros(N,N);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    for i = 1:N
        waitbar(i/N,f)
        for j = 1:N
            x0 = [X(i),Y(j)];
            [~,x] = ode45(@nlsys,[0,T],x0,opts);
            D = zeros(Npf,2);
            for k = 1:Npf
                D(k,:) = [k,norm(pf(k,:)-x(end,:))];
            end   
            idx = D(:,2) == min(D(:,2));
            B(j,i) = D(idx,1);
        end
    end
    delete(f)
end
%
%
%
function dxdt = nlsys(t,x)
    dxdt = [x(2);-100*x(2)-8*(x(1)^3)+2*x(1)];
end