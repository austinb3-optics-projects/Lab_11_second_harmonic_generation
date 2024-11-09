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

zspan = [0, 1e-2];
initial = [A1, 0];

% define a function that will represent the system of equations
function dAdz = system(z,A,eta)

    dAdz = zeros(2,1);
    % specify the diff eqs in a column vector
    dAdz(1) = 1i*eta*A(2)*conj(A(1));
    dAdz(2) = 1i*eta*A(1)^2;


end
% since more than just the z,A are needed to define the diff eqs, the ode
% function must be a function of the function of diff eqs such that we can
% input what that constant must be. In this case, eta.
odefun = @(z,A) system(z,A,eta);
options = odeset('RelTol',1e-7,'AbsTol',1e-8,'MaxStep',0.01);


[z,A] = ode45(odefun,zspan,initial,options);

n2 = n1;        % true in the phase matched case.

I_1 = 2*eps0*n1*c*(A(:,1).*conj(A(:,1)));
I_2 = 2*eps0*n2*c*(A(:,2).*conj(A(:,2)));


figure;
plot(z,I_1,'Linewidth',1.5); hold on;
plot(z,I_2,'LineWidth',1.5); hold off; grid on;
title('Fundamental and SH intensity');
ylabel('$I_{1,2}$','Interpreter','latex');
xlabel('$z$','Interpreter','latex');
set(gca,'FontSize',15);
legend('I_{1}','I_{2}');
