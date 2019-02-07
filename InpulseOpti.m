clear all
close all
%% Initalise System's Parameters
k = 5;   % spring constant
b = 0.5; % damping coefficient
% System Matrices
A = [0,1;-k,-b];
B = [0;1];
% Initial State
x0 = [-0.8,0];

%% Optimization setup
options = optimoptions('fmincon',...
    'Display','iter',...
    'Algorithm','sqp',...
    'MaxFunctionEvaluations',1e3);%,...
%    'PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt});
%
opts = odeset('RelTol',1e-3,'AbsTol',1e-4);
%% 1st Trajectory (x0 --> x1*) 
[t1,x1] = ode45(@nlsys_noisy,[0,5],x0,opts);
xi = x1(end,:)';
xj = [.5;0];

%% Optimise Impulse #1 (x1* --> x2*)
% Run Optimizer
tic
var = fmincon(@(var) cost(var,A,B,xi,xj),[1,1],[],[],[],[],...
    [0.1,-1000],[100,1000],@(var) basincon(var,A,B,xi,xj),options);
toc
td = var(1);
gamma = var(2);
% Check Result
xj_a = expm(td*A)*(xi+B*gamma);
% Compute Trajectory
t2 = 0:0.001:td;
x2 = zeros(length(t2),2);
%
for i = 1:length(t2)
   x2(i,:) = (expm(t2(i)*A)*(xi+B*gamma))';
end

%% 2nd Trajectory (x1d -->x2*)
% New Initial Condition
x0 = x2(end,:);
%
[t3,x3] = ode45(@nlsys_noisy,[0,5],x0,opts);
%
xi = x3(end,:)';
xj = [-.5;0];

%% Optimise Impulse #2 (x2* --> x1*)
% Run Optimizer
tic
var2 = fmincon(@(var2) cost(var2,A,B,xi,xj),[1,1],[],[],[],[],...
    [0.1,-1000],[100,1000],@(var) basincon(var,A,B,xi,xj),options);
toc
gamma = var2(2);
td = var2(1);
% Check Result
xj_a = expm(td*A)*(xi+B*gamma);
% Compute Trajectory
t4 = 0:0.001:td;
x4 = zeros(length(t4),2);
for i = 1:length(t4)
   x4(i,:) = (expm(t4(i)*A)*(xi+B*gamma))';
end

%% 3rd Trajectory (x2d --> x1*)
x0 = x4(end,:);
[t5,x5] = ode45(@nlsys_noisy,[0,5],x0,opts);
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute error norms
xj = [0.5;0];
xi = [-0.5;0];
e1 = xj - expm(var(1)*A)*(xi+B*var(2));
ne1 = norm(e1)
xj = [-0.5;0];
xi = [0.5;0];
e2 = xj - expm(var2(1)*A)*(xi+B*var2(2));
ne2 = norm(e2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results
% Append state & time vectors
t = [t1;t2'+t1(end);t3+t1(end)+t2(end);t4'+t3(end)+t1(end)+t2(end);t5+t4(end)+t3(end)+t1(end)+t2(end)];
x = [x1;x2;x3;x4;x5];
% Plots
% Time Evolution
figure()
subplot(211)
hold on
plot([t1(1),t1(end)+t2(end)+t3(end)+t4(end)+t5(end)],[0.5,0.5],'--g','LineWidth',1)
plot([t1(1),t1(end)+t2(end)+t3(end)+t4(end)+t5(end)],[-0.5,-0.5],'--g','LineWidth',1)
plot(t1,x1(:,1),'k', 'LineWidth',1.5)
plot(t1(end)+t2,x2(:,1),'-.k', 'LineWidth',1.5)
plot(t1(end)+t2(end)+t3,x3(:,1),'k', 'LineWidth',1.5)
plot(t1(end)+t2(end)+t3(end)+t4,x4(:,1),'-.k', 'LineWidth',1.5)
plot(t1(end)+t2(end)+t3(end)+t4(end)+t5,x5(:,1),'k', 'LineWidth',1.5)
plot(t1(end)+t2(end)+t3(end),x3(end,1),'--ob',...
    'MarkerEdgeColor',[0.,0.,0.]','MarkerFaceColor',[0.5,0.5,0.5]','MarkerSize',4.5,'LineWidth',1.5)
plot(t1(end),x1(end,1),'--ob',...
    'MarkerEdgeColor',[0.,0.,0.]','MarkerFaceColor',[0.5,0.5,0.5]','MarkerSize',4.5,'LineWidth',1.5)
hold off
xlabel('$t$ [$s$]')
ylabel('$\xi(t)$')
box on
subplot(212)
hold on
plot([t1(1),t1(end)+t2(end)+t3(end)+t4(end)+t5(end)],[0.,0.],'--g','LineWidth',1)
plot(t1,x1(:,2),'k', 'LineWidth',1.5)
plot(t1(end)+t2,x2(:,2),'-.k', 'LineWidth',1.5)
%plot([t1(end),t1(end)],[x2(end,2),x3(1,2)])
plot(t1(end)+t2(end)+t3,x3(:,2),'k', 'LineWidth',1.5)
plot(t1(end)+t2(end)+t3(end)+t4,x4(:,2),'-.k', 'LineWidth',1.5)
plot(t1(end)+t2(end)+t3(end)+t4(end)+t5,x5(:,2),'k', 'LineWidth',1.5)
plot([t1(end),t1(end)],[x1(end,2),x2(1,2)],':ob',...
    'MarkerEdgeColor',[0.,0.,0.]','MarkerFaceColor',[0.5,0.5,0.5]','MarkerSize',4.5,'LineWidth',1.5)
plot([t1(end)+t2(end)+t3(end),t1(end)+t2(end)+t3(end)],[x3(end,2),x4(1,2)],':ob',...
    'MarkerEdgeColor',[0.,0.,0.]','MarkerFaceColor',[0.5,0.5,0.5]','MarkerSize',4.5,'LineWidth',1.5)
hold off
xlabel('$t$ [$s$]')
ylabel('$\dot{xi}(t)$')
box on
%
% State Space Trajectory
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
xlabel('$\xi(t)$')
ylabel('$\dot{xi}(t)$')
box on
%
basincon(var2,A,B,xi,xj)
xff = limsol(var2,A,B,xi);
%
% Comoute Control efforts
% energy shaping
u_es = [-8*(x1(:,1).^3)+(2+5)*x1(:,1);...
       zeros(length(t2),1);...
       -8*(x3(:,1).^3)+(2+5)*x3(:,1);...
       zeros(length(t4),1);...
       -8*(x5(:,1).^3)+(2+5)*x5(:,1)];

u_di = [-4.5.*x1(:,2);...
       zeros(length(t2),1);...
       -4.5.*x3(:,2);...
       zeros(length(t4),1);...
       -4.5.*x5(:,2)];
figure()
box on
plot(t,u_es,'k','LineWidth',2)
hold on
plot(t,u_di,':k','LineWidth',2)
xlabel('t')
ylabel('\beta(x(t)), v(t)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%
% Cost Function
function IDX = cost(var,A,B,xi,xj)
    e = xj - expm(var(1)*A)*(xi+B*var(2));
    IDX = 0.001*var(2)^2 + norm(e)^2;
    %IDX = IDX + ((xi+B*var(2))')*( ((A+A')^-1)*(expm(var(1)*(A+A'))-eye(2)))*(xi+B*var(2));
end
% Nonlinear Constraint (Appartenenza al basin of attraction)
function [c,ceq] = basincon(var,A,B,xi,xj)
    c = norm(xj-limsol(var,A,B,xi)') - 1e-5;
    ceq = [];
end
% Compute inf. time limit of the trajectory (lim_{t->inf} Phi(t,x0,u))
function xl = limsol(var,A,B,xi)
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    xtd = expm(var(1)*A)*(xi+B*var(2));
    [~,x] = ode45(@nlsys,[0,1000],xtd');
    xl = x(end,:);
end
% EB-PBC Controlled System (Nonlinear)
function dxdt = nlsys(t,x)
    dxdt = [x(2);-5*x(2)-8*(x(1)^3)+2*x(1)];
end
% EB-PBC Controlled System + NOISE ON STATE MEASURENTS
function dxdt = nlsys_noisy(t,x)
    std = 0.0;
    x = x + std*randn(2,1);
    dxdt = [x(2);-5*x(2)-8*(x(1)^3)+(2)*x(1)];
end

