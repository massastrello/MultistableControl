clear all
close all
%
k = 5;
b = .5;
%
A = [0,1;-k,-b];
B = [0;1];
%
x0 = [-0.8,0];
%
[t1,x1] = ode45(@nlsys,[0,15],x0);
xi = x1(end,:)';
xj = [.5;0];

options = optimoptions('fmincon',...
    'Display','iter',...
    'Algorithm','sqp',...
    'MaxFunctionEvaluations',1e3);%,...
%    'PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt});

%
tic
var = fmincon(@(var) cost(var,A,B,xi,xj),[1,1],[],[],[],[],...
    [0.1,-1000],[100,1000],@(var) basincon(var,A,B,xi,xj),options);
toc
% Vallidation
td = var(1);
gamma = var(2);
xj_a = expm(td*A)*(xi+B*gamma);


t2 = 0:0.001:td;
x2 = zeros(length(t2),2);
%
for i = 1:length(t2)
   x2(i,:) = (expm(t2(i)*A)*(xi+B*gamma))';
end
%
%
x0 = x2(end,:);
%
[t3,x3] = ode45(@nlsys,[0,15],x0);
%
%
%
xi = x3(end,:)';%[-.5;0];
xj = [-.5;0];
tic
var2 = fmincon(@(var2) cost(var2,A,B,xi,xj),[1,1],[],[],[],[],...
    [0.1,-1000],[100,1000],@(var) basincon(var,A,B,xi,xj),options);
toc
gamma = var2(2);
td = var2(1);
xj_a = expm(td*A)*(xi+B*gamma);
t4 = 0:0.001:td;
x4 = zeros(length(t4),2);
for i = 1:length(t4)
   x4(i,:) = (expm(t4(i)*A)*(xi+B*gamma))';
end
%
%
x0 = x4(end,:);
[t5,x5] = ode45(@nlsys,[0,15],x0);
%
%
t = [t1;t2'+t1(end);t3+t1(end)+t2(end);t4'+t3(end)+t1(end)+t2(end);t5+t4(end)+t3(end)+t1(end)+t2(end)];
x = [x1;x2;x3;x4;x5];
%
figure()
plot(t,x)
%
%
figure(2)
hold on 
plot(x1(:,1),x1(:,2),'k', 'LineWidth',2)
plot([x1(end,1),x2(1,1)],[x1(end,2),x2(1,2)],':b', 'LineWidth',2)
plot(x2(:,1),x2(:,2),'-.k', 'LineWidth',2)
plot([x2(end,1),x3(1,1)],[x2(end,2),x3(1,2)],':b', 'LineWidth',2)
plot(x3(:,1),x3(:,2),'k', 'LineWidth',2)
plot([x3(end,1),x4(1,1)],[x3(end,2),x4(1,2)],':b', 'LineWidth',2)
plot(x4(:,1),x4(:,2),'-.k', 'LineWidth',2)
plot([x4(end,1),x5(1,1)],[x4(end,2),x5(1,2)],':b', 'LineWidth',2)
plot(x5(:,1),x5(:,2),'k', 'LineWidth',2)
plot([x1(end,1),x2(1,1)],[x1(end,2),x2(1,2)],'ok', 'MarkerFaceColor',[0.5,0.5,0.5]','LineWidth',1)
plot([x3(end,1),x4(1,1)],[x3(end,2),x4(1,2)],'ok', 'MarkerFaceColor',[0.5,0.5,0.5]','LineWidth',1)
hold off

basincon(var2,A,B,xi,xj)
xff = limsol(var2,A,B,xi);
%
function IDX = cost(var,A,B,xi,xj)
    e = xj - expm(var(1)*A)*(xi+B*var(2));
    IDX = 0.001*var(2)^2 + norm(e)^2;
    %IDX = IDX + ((xi+B*var(2))')*( ((A+A')^-1)*(expm(var(1)*(A+A'))-eye(2)))*(xi+B*var(2));
end
%
function [c,ceq] = basincon(var,A,B,xi,xj)
    c = norm(xj-limsol(var,A,B,xi)') - 1e-5;
    ceq = [];
end
%
function xl = limsol(var,A,B,xi)
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    xtd = expm(var(1)*A)*(xi+B*var(2));
    [~,x] = ode45(@nlsys,[0,1000],xtd');
    xl = x(end,:);
end

function dxdt = nlsys(t,x)
    dxdt = [x(2);-5*x(2)-8*(x(1)^3)+2*x(1)];
end

