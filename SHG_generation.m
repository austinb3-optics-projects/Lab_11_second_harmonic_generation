clear; close all;

%% initialization


deff = 1e-12;   % [m/V]: effective nonlinear constant.
lam1 = 1e-6;    % wavelength of the fundamental
n1 = 1.6;       % index of the fundamental.
c = 3e8;

%% 1a)
k1 = 2*pi/lam1;
w1 = k1*c;

eta = 2*k1*deff;
%% 1b)

I1 = 1e13;  % peak intensity of a fundamental gaussian beam at r = 0, z=0
eps0 = 8.8541e-12;      % permittivity of free space.
w0 = 10e-6;             % beam waist

A1 = sqrt(I1/2/eps0/n1/c);
P1 = pi*eps0*n1*c*A1^2 * w0^2;

%% Part B

I1 = 1e13; 
L = 1e-2;

zspan = [0, L];
initial = [A1, 0];
dk = 0;     % perfectly phase matched
% define a function that will represent the system of equations
function dAdz = coupled_eqs(z,A,eta,dk)

    dAdz = zeros(2,1);
    % specify the diff eqs in a column vector
    dAdz(1) = 1i*eta*A(2)*conj(A(1))*exp(-1i*dk*z);
    dAdz(2) = 1i*eta*A(1)^2*exp(1i*dk*z);


end
% since more than just the z,A are needed to define the diff eqs, the ode
% function must be a function of the function of diff eqs such that we can
% input what that constant must be. In this case, eta.
odefun = @(z,A) coupled_eqs(z,A,eta,dk);
options = odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',0.01);


[z,A] = ode45(odefun,zspan,initial,options);

n2 = n1;        % true in the phase matched case.

I_1 = 2*eps0*n1*c*(A(:,1).*conj(A(:,1)));
I_2 = 2*eps0*n2*c*(A(:,2).*conj(A(:,2)));

% numeric
figure(1);
plot(z,I_1,'Linewidth',1.5); hold on;
plot(z,I_2,'LineWidth',1.5); hold off; grid on;
title('Fundamental and SH intensity');
ylabel('$I_{1,2}$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
set(gca,'FontSize',15);
legend('I_{1}','I_{2}');

%% part 2b analytical solutions 

l = sqrt(2*n1^2 *n2 * eps0 * c^3)/(2*w1*deff*sqrt(I1));

A1_an = A1.*sech(z./l);
A2_an = A1.*tanh(z./l);

I_1_an = 2*eps0*n1*c*(A1_an.*conj(A1_an));
I_2_an = 2*eps0*n1*c*(A2_an.*conj(A2_an));

% analytic
figure(2);
plot(z,I_1_an,'Linewidth',1.5); hold on;
plot(z,I_2_an,'LineWidth',1.5); hold off; grid on;
title('Fundamental and SH intensity (analytic solution)');
ylabel('$I_{1,2}$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
set(gca,'FontSize',15);
legend('I_{1}','I_{2}');

% analytic vs numeric
figure(3);
plot(z,I_1,'Linewidth',1.5); hold on;
plot(z,I_2,'LineWidth',1.5); grid on;
plot(z,I_1_an,'Linewidth',1.5); 
plot(z,I_2_an,'LineWidth',1.5); hold off;
title('Analytic vs Numeric');
ylabel('$I_{1,2}$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
set(gca,'FontSize',15);
legend('I_{1,num}','I_{2,num}','I_{1,an}','I_{2,an}');


%% part 2c phase mismatch

I1 = 0.1e13;    % incident intensity
A1 = sqrt(I1/2/eps0/n1/c);
initial = [A1, 0];
% phase mismatch choices 
dn = [1e-4,3e-4];
dk = 2*w1*dn/c;
Lcoh = pi./abs(dk);
Period = 2.*Lcoh;


odefun1 = @(z,A) coupled_eqs(z,A,eta,dk(1));
opts = odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',0.00001);

[zn1,An1] = ode45(odefun1,zspan,initial,opts);

In1 = 2*eps0*n1*c*(An1(:,1).*conj(An1(:,1)));
In12 = 2*eps0*n2*c*(An1(:,2).*conj(An1(:,2)));


odefun2 = @(z,A) coupled_eqs(z,A,eta,dk(2));

[zn2,An2] = ode45(odefun2,zspan,initial,opts);


In21 = 2*eps0*n1*c*(An2(:,1).*conj(An2(:,1)));
In2 = 2*eps0*n2*c*(An2(:,2).*conj(An2(:,2)));

% need to plot separately because the zni parameter is different between
% each solution

% dn = 1e-4;
figure(4);
plot(zn1,In12/I1,'LineWidth',1.5); grid on;
title('Phase mismatched case');
ylabel('$I_{2}(z)/I_{1}(0)$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
legend('\Deltan = 10^{-4}');
set(gca,'FontSize',15);

% dn = 3e-4
figure(5);
plot(zn2,In2/I1,'LineWidth',1.5); grid on;
title('Phase mismatched case');
ylabel('$I_{2}(z)/I_{1}(0)$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
legend('\Deltan = 3\times 10^{-4}');
set(gca,'FontSize',15);

%% part 2d Undepleted pump beam approximation

% this section deals with an approximation in which A_1(z) = A_1(0) 

z = linspace(0,L,300);

intfun = @(z) 1i*eta*A1^2 *exp(1i*dk(2)*z);
A2z = zeros(1,length(z));
for i = 1:length(z)
    A2z(i) = integral(intfun,0,z(i),'WayPoints', [pi/2, pi, 3*pi/2, 2*pi]);
end

I2z = 2*eps0*n2*c*(A2z.*conj(A2z));

figure(6);
plot(z,I2z./I1,'LineWidth',1.5); grid on;
title('undepleted pump appox');
ylabel('$I_{2}(z)/I_{1}(0)$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
legend('\Deltan = 3\times 10^{-4}');
set(gca,'FontSize',15);

%% Part C Gaussian beam propagation

% expressions : b = k1*w0^2
% wz = w0*sqrt(b^2 + 4z^2)/b
%phi = arctan(2z/b)

w0 = 10e-6;     % waist size
b = k1*w0^2;    % rayleigh range.

z = linspace(-3*b,3*b,400);

wz = w0.*sqrt(b^2 + 4.*z.^2)./b;
phi = atan(2.*z./b);


figure(7);
plot(z,wz./w0,'LineWidth',1.5); grid on;
title('Beam size');
ylabel('$w(z)/w_0$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
set(gca,'FontSize',15);


figure(8);
plot(z,phi./z,'LineWidth',1.5); grid on;
title('Guoy phase');
ylabel('$\Phi(z)/\pi$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
set(gca,'FontSize',15);

%% Section II 
% a and b
z0 = -L/2;
z_last = L/2;       % suggests that the gaussian beam focuses at the center of the crystal.

% xi = L/b; 
% sigma = b*dk/2;
% zeta = 2.*z./b;

% SH waist size is W2 = W1./sqrt(2)
% b = k1*w1^2
% b = k2*w2^2
% a factor of 2 comes from k2 = 2k1 and a factor of 1/2 comes from w2^2 =
% w1^2/2

xi = 2.84;
sigma = linspace(-3,3,400);



h = zeros(1,length(sigma));
for i = 1:length(sigma)
    intfun = @(zeta) exp(1i.*sigma(i).*zeta)./(1 + 1i*zeta);
    h(i) = (1/4.27/xi).*abs(integral(intfun,-xi,xi)).^2;
end

figure(9);
plot(sigma,h,'LineWidth',1.5); grid on;
title ('Geometric factor');
ylabel('$h(\sigma,\xi)$','Interpreter','latex');
xlabel('$\sigma$','Interpreter','latex');
set(gca,'FontSize',15);
%% Section II part c

xi = linspace(1,6,400);
sigma = 0.56;

intfun = @(zeta) exp(1i.*sigma.*zeta)./(1 + 1i*zeta);

h = zeros(1,length(xi));
for i = 1:length(xi)
    h(i) = (1/4.27/xi(i)).*abs(integral(intfun,-xi(i),xi(i))).^2;
end

figure(9);
plot(xi,h,'LineWidth',1.5); grid on;
title ('Geometric factor');
ylabel('$h(\sigma,\xi)$','Interpreter','latex');
xlabel('$\xi$','Interpreter','latex');
set(gca,'FontSize',15);

%% section II part d

xi = linspace(1,4,400);
sigma = linspace(0,2,400);

h = zeros(length(sigma),length(xi));
for i = 1:length(xi)
    for j = 1:length(sigma)
        intfun = @(zeta) exp(1i.*sigma(j).*zeta)./(1 + 1i*zeta);
        h(j,i) = (1/4.27/xi(i)).*abs(integral(intfun,-xi(i),xi(i))).^2;
    end
end

figure(10);
imagesc([xi(1),xi(end)],[sigma(1),sigma(end)],h);
colormap jet;
colorbar;
axis xy;
title ('Geometric factor');
ylabel('$\xi$','Interpreter','latex');
xlabel('$\sigma$','Interpreter','latex');
set(gca,'FontSize',15);

%% section II part e

% expressions : b = k1*w0^2
% wz = w0*sqrt(b^2 + 4z^2)/b
%phi = arctan(2z/b)
L = 3e-3;

xi = L/b;

z = linspace(-L/2,L/2,201);

sigma = linspace(-3,3,400);



h = zeros(1,length(sigma));
for i = 1:length(sigma)
    intfun = @(zeta) exp(1i.*sigma(i).*zeta)./(1 + 1i*zeta);
    h(i) = (1/4.27/xi).*abs(integral(intfun,-xi,xi)).^2;
end
[hmax, idx] = max(h);
sigma_max = sigma(idx);
dk = 2*sigma_max/b;
dn = dk/2/k1;
fprintf('The required sigma = %f and the required dn = %f\n',sigma_max,dn);
% % % figure;
% % % plot(sigma,h,'LineWidth',1.5); grid on;
% % % title ('Geometric factor');
% % % ylabel('$h(\sigma,\xi)$','Interpreter','latex');
% % % xlabel('$\sigma$','Interpreter','latex');
% % % set(gca,'FontSize',15);

%% Section II part f

deff = 15e-12;

P1 = 10;

% P2 = KLk1hP1^2

K = 2.14*w1.^2*deff.^2./eps0./n1.^2./n2./c.^3./pi;
ENL = K*L*hmax*k1;
P2 = ENL.*P1.^2;

fprintf('the conversion efficiency is %f\n',ENL*100);
fprintf('The SH output power is P_2 = %f W\n',P2);