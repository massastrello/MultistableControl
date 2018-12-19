clear all
close all
%
k = 5;
b = 0;
%
A = [0,1;-k,-b];
B = [0;1];
%
xi = [-.5;0];
xj = [.5;0];
%
l = [0,-3];
u = [3,3];
N = 100;

b = linspace(0,3,6);
for i = 1:6
    A = [0,1;-k,-b(i)];
    [td,gamma,J] = CostMap(l,u,N,A,B,xi,xj);
    Ji(:,:,i) = J;
    [Xm,Ym] = meshgrid(td,gamma);
    figure()
    surf(Xm,Ym,J,'EdgeColor','None')
    view(2)
end

[Xm,Ym] = meshgrid(td,gamma);
%Bm = griddata(X,Y,B,Xm,Ym,'natural');

figure()
surf(Xm,Ym,J,'EdgeColor','None')
view(2)

function [td,gamma,J] = CostMap(l,u,N,A,B,xi,xj)
    %f = waitbar(0,'Please wait...');
    td = linspace(l(1),u(1),N);
    gamma = linspace(l(2),u(2),N);
    %
    J = zeros(N,N);
    parfor i = 1:N
        %waitbar(i/N,f)
        for j = 1:N
            var = [td(i),gamma(j)];
            J(j,i) = cost(var,A,B,xi,xj)
        end
    end
    %delete(f)
end

function IDX = cost(var,A,B,xi,xj)
    e = xj - expm(var(1)*A)*(xi+B*var(2));
    IDX = norm(e)^2 ;%+ var(2)^2;
    %IDX = IDX + ((xi+B*var(2))')*( ((A+A')^-1)*(expm(var(1)*(A+A'))-eye(2)))*(xi+B*var(2));
end