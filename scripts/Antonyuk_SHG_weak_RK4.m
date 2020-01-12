%Generacja drugiej harmonicznej o niskiej intesywnosci 
%na podstawie parametrów z 
%Boris P. Antonyuk, Light-Driven Alignment, Springer, str. 118–123, 2009.
%
% Autor: Sylwia Majchrowska
% 7 kwietnia 2017r.

clear all;
clc;
format long e

M = 500;
N = 500;
L = 5;
time = 100;

s = linspace(0, L, N);
tau = linspace(0, time, M);
dS = s(2) - s(1);
dtau = tau(2) - tau(1);

E0 = zeros(M, N);
E2 = zeros(M, N);
Psi0 = zeros(M, N);
Psi2 = zeros(M, N);

% warunki poczatkowe
E0(1,:) = 1e-15;  % E0(S,0)
Psi0(1,:) = 0;  % Psi0(S,0)
E2(:,1) = 1e-15;  % E2(0,tau)
Psi2(:,1) = 0;  % Psi2(0,tau)

%Psi2(tau>25,1) = pi;  % zmiana fazy dla tau = 25

rhs1 = @(E0,Psi0,Psi2)    E0 * sin(Psi0 - Psi2);
rhs2 = @(E0,E2,Psi0,Psi2) -E0/E2 * cos(Psi0 - Psi2);
rhs3 = @(E0,E2,Psi0,Psi2) -E2^2 *(E0 - cos(Psi0 - Psi2));
rhs4 = @(E0,E2,Psi0,Psi2) -E2^2/E0 * sin(Psi0 - Psi2);

fprintf(1, '\nStart...       ');
tic

for j = 1:N
    for i = 1:M
        if(j ~= N)
            kE2_1   = dS * rhs1(E0(i,j),                   Psi0(i,j),Psi2(i,j));
            kPsi2_1 = dS * rhs2(E0(i,j),E2(i,j),           Psi0(i,j),Psi2(i,j));
            kE2_2   = dS * rhs1(E0(i,j),                   Psi0(i,j),Psi2(i,j) + .5*kPsi2_1);
            kPsi2_2 = dS * rhs2(E0(i,j),E2(i,j) + .5*kE2_1,Psi0(i,j),Psi2(i,j) + .5*kPsi2_1);
            kE2_3   = dS * rhs1(E0(i,j),                   Psi0(i,j),Psi2(i,j) + .5*kPsi2_2);
            kPsi2_3 = dS * rhs2(E0(i,j),E2(i,j) + .5*kE2_2,Psi0(i,j),Psi2(i,j) + .5*kPsi2_2);
            kE2_4   = dS * rhs1(E0(i,j),                   Psi0(i,j),Psi2(i,j) + kPsi2_3);
            kPsi2_4 = dS * rhs2(E0(i,j),E2(i,j) + kE2_3,   Psi0(i,j),Psi2(i,j) + kPsi2_3);
            E2(i,j+1)   = E2(i,j)   + (kE2_1   + 2*kE2_2   + 2*kE2_3   + kE2_4   ) / 6;
            Psi2(i,j+1) = Psi2(i,j) + (kPsi2_1 + 2*kPsi2_2 + 2*kPsi2_3 + kPsi2_4 ) / 6;
        end
        if(i ~= M)
            kE0_1   = dtau * rhs3(E0(i,j),           E2(i,j),Psi0(i,j),             Psi2(i,j));
            kPsi0_1 = dtau * rhs4(E0(i,j),           E2(i,j),Psi0(i,j),             Psi2(i,j));
            kE0_2   = dtau * rhs3(E0(i,j) + .5*kE0_1,E2(i,j),Psi0(i,j) + .5*kPsi0_1,Psi2(i,j));
            kPsi0_2 = dtau * rhs4(E0(i,j) + .5*kE0_1,E2(i,j),Psi0(i,j) + .5*kPsi0_1,Psi2(i,j));
            kE0_3   = dtau * rhs3(E0(i,j) + .5*kE0_2,E2(i,j),Psi0(i,j) + .5*kPsi0_2,Psi2(i,j));
            kPsi0_3 = dtau * rhs4(E0(i,j) + .5*kE0_2,E2(i,j),Psi0(i,j) + .5*kPsi0_2,Psi2(i,j));
            kE0_4   = dtau * rhs3(E0(i,j) + kPsi0_3, E2(i,j),Psi0(i,j) + kPsi0_3,   Psi2(i,j));
            kPsi0_4 = dtau * rhs4(E0(i,j) + kE0_3,   E2(i,j),Psi0(i,j) + kPsi0_3,   Psi2(i,j));
            E0(i+1,j)   = E0(i,j)   + (kE0_1   +  2*kE0_2   + 2*kE0_3   + kE0_4   ) / 6;
            Psi0(i+1,j) = Psi0(i,j) + (kPsi0_1 + 2*kPsi0_2 + 2*kPsi0_3 + kPsi0_4 ) / 6;
        end
        fprintf(1, '\b\b\b\b\b\b\b%06.2f%%', (i + M*(j-1)) * 100.0 /M/N );
    end
end

tx = toc;

fprintf(1, '\n\nCzas trwania symulacji (s) = ');
fprintf(1, '%5.2f%', tx );
fprintf(1, '\n\n');

%--------------------------------------------------------------------------
%Wykresy
%--------------------------------------------------------------------------

figure(1)
set(gca, 'fontsize', 14)
plot(s, E0(M,:), 'g', s, E2(M,:), 'r')  % dla tau = 100
grid on
legend('E_0', 'E_2', 'Location', 'southeast')
xlabel('S'); ylabel('E_i');
% print('-f1','E_week_RK4','-dpng')

figure(2)
set(gca, 'fontsize', 14)
plot(s, Psi0(M,:), 'g', s, Psi2(M,:), 'r')
grid on
legend('\Psi_0', '\Psi_2', 'Location', 'northeast')
xlabel('S'); ylabel('\Psi_i');
% print('-f2','Psi_week_RK4','-dpng')

figure(3)
set(gca, 'fontsize', 14)
set(gcf,'renderer','zbuffer');
mesh(s, tau, E2)
set(gca, 'xlim', [0 L], 'ylim', [0 time])
xlabel('S'); ylabel('\tau'); zlabel('E_2');
% print('-f3','E2_week_mesh_RK4', '-dpng')

figure(4)
set(gca, 'fontsize', 14)
set(gcf,'renderer','zbuffer');
mesh(s, tau, E2.^2)
set(gca, 'xlim', [0 L], 'ylim', [0 time])
xlabel('S'); ylabel('\tau'); zlabel('E_2');
% print('-f4','E2_week_mesh_RK4', '-dpng')

% close all
% clc