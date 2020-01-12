%Generacja drugiej harmonicznej o wysokiej intesywnosci 
%na podstawie parametrów z 
%Boris P. Antonyuk, Light-Driven Alignment, str. 123–128, Springer 2009.
%
% Autor: Sylwia Majchrowska
% 7 kwietnia 2017r.

clear all;
clc;
format long e

M = 2^17;
N = 2^10;
L = 10;
time = 5000;

s = linspace(0, L, N);
tau = linspace(0, time, M);
dS = s(2) - s(1);
dtau = tau(2) - tau(1);

K = 500;
E0 = zeros(M, N);
E1 = zeros(M, N);
E2 = zeros(M, N);
Psi0 = zeros(M, N);
Psi1 = zeros(M, N);
Psi2 = zeros(M, N);

% warunki poczatkowe
E0(1,:) = 1e-15;  % E0(S,0)
Psi0(1,:) = 0;  % Psi0(S,0)
E2(:,1) = 0.05;  % E2(0,tau)
Psi2(:,1) = 0;  % Psi2(0,tau)
E1(:,1) = 1;  % E1(0,tau)
Psi1(:,1) = 0;  % Psi1(0,tau)

rhs1 = @(E0,E1,Psi0,Psi1,Psi2)    E0 * E1^2 * sin(Psi0 + 2 * Psi1 - Psi2);
rhs2 = @(E0,E1,E2,Psi0,Psi1,Psi2) -E0 * E1^2/E2 * cos(Psi0 + 2 * Psi1 - Psi2);
rhs3 = @(E0,E1,E2,Psi0,Psi1,Psi2) -E0 * E1 * E2 * sin(Psi0 + 2 * Psi1 - Psi2);
rhs4 = @(E0,E2,Psi0,Psi1,Psi2)    -E0 * E2 * cos(Psi0 + 2 * Psi1 - Psi2);
rhs5 = @(E0,E2,Psi0,Psi1,Psi2)    -E2^2 *(E0 - cos(Psi0 + 2 * Psi1 - Psi2));
rhs6 = @(E0,E2,Psi0,Psi1,Psi2)    -E2^2/E0 * sin(Psi0 + 2 * Psi1 - Psi2);

fprintf(1, '\nStart...       ');
tic
for j = 1:N
    for i = 1:M
        if(j ~= N)
            kE2_1   = dS * rhs1(E0(i,j),E1(i,j),                   Psi0(i,j),Psi1(i,j),Psi2(i,j));
            kPsi2_1 = dS * rhs2(E0(i,j),E1(i,j),E2(i,j),           Psi0(i,j),Psi1(i,j),Psi2(i,j));
            kE1_1   = dS * rhs3(E0(i,j),E1(i,j),E2(i,j),           Psi0(i,j),Psi1(i,j),Psi2(i,j));
            kPsi1_1 = dS * rhs4(E0(i,j),        E2(i,j),           Psi0(i,j),Psi1(i,j),Psi2(i,j));
            
            kE2_2   = dS * rhs1(E0(i,j),E1(i,j) + .5*kE1_1,                   Psi0(i,j),Psi1(i,j) + .5*kPsi1_1,Psi2(i,j) + .5*kPsi2_1);
            kPsi2_2 = dS * rhs2(E0(i,j),E1(i,j) + .5*kE1_1,E2(i,j) + .5*kE2_1,Psi0(i,j),Psi1(i,j) + .5*kPsi1_1,Psi2(i,j) + .5*kPsi2_1);
            kE1_2   = dS * rhs3(E0(i,j),E1(i,j) + .5*kE1_1,E2(i,j) + .5*kE2_1,Psi0(i,j),Psi1(i,j) + .5*kPsi1_1,Psi2(i,j) + .5*kPsi2_1);
            kPsi1_2 = dS * rhs4(E0(i,j),                   E2(i,j) + .5*kE2_1,Psi0(i,j),Psi1(i,j) + .5*kPsi1_1,Psi2(i,j) + .5*kPsi2_1);
            
            kE2_3   = dS * rhs1(E0(i,j),E1(i,j) + .5*kE1_2,                   Psi0(i,j),Psi1(i,j) + .5*kPsi1_2,Psi2(i,j) + .5*kPsi2_2);
            kPsi2_3 = dS * rhs2(E0(i,j),E1(i,j) + .5*kE1_2,E2(i,j) + .5*kE2_2,Psi0(i,j),Psi1(i,j) + .5*kPsi1_2,Psi2(i,j) + .5*kPsi2_2);
            kE1_3   = dS * rhs3(E0(i,j),E1(i,j) + .5*kE1_2,E2(i,j) + .5*kE2_2,Psi0(i,j),Psi1(i,j) + .5*kPsi1_2,Psi2(i,j) + .5*kPsi2_2);
            kPsi1_3 = dS * rhs4(E0(i,j),                   E2(i,j) + .5*kE2_2,Psi0(i,j),Psi1(i,j) + .5*kPsi1_2,Psi2(i,j) + .5*kPsi2_2);           
            
            kE2_4   = dS * rhs1(E0(i,j),E1(i,j) + kE1_3,                Psi0(i,j),Psi1(i,j) + kPsi1_3,Psi2(i,j) + kPsi2_3);
            kPsi2_4 = dS * rhs2(E0(i,j),E1(i,j) + kE1_3,E2(i,j) + kE2_3,Psi0(i,j),Psi1(i,j) + kPsi1_3,Psi2(i,j) + kPsi2_3);
            kE1_4   = dS * rhs3(E0(i,j),E1(i,j) + kE1_3,E2(i,j) + kE2_3,Psi0(i,j),Psi1(i,j) + kPsi1_3,Psi2(i,j) + kPsi2_3);
            kPsi1_4 = dS * rhs4(E0(i,j),                E2(i,j) + kE2_3,Psi0(i,j),Psi1(i,j) + kPsi1_3,Psi2(i,j) + kPsi2_3);
            
            E2(i,j+1)   = E2(i,j)   + (kE2_1   + 2*kE2_2   + 2*kE2_3   + kE2_4   ) / 6;
            Psi2(i,j+1) = Psi2(i,j) + (kPsi2_1 + 2*kPsi2_2 + 2*kPsi2_3 + kPsi2_4 ) / 6;
            E1(i,j+1)   = E1(i,j)   + (kE1_1   + 2*kE1_2   + 2*kE1_3   + kE1_4   ) / 6;
            Psi1(i,j+1) = Psi1(i,j) + (kPsi1_1 + 2*kPsi1_2 + 2*kPsi1_3 + kPsi1_4 ) / 6;
            
        end
        if(i ~= M)
            kE0_1   = dtau * rhs5(E0(i,j),           E2(i,j),Psi0(i,j),             Psi1(i,j),Psi2(i,j));
            kPsi0_1 = dtau * rhs6(E0(i,j),           E2(i,j),Psi0(i,j),             Psi1(i,j),Psi2(i,j));
            kE0_2   = dtau * rhs5(E0(i,j) + .5*kE0_1,E2(i,j),Psi0(i,j) + .5*kPsi0_1,Psi1(i,j),Psi2(i,j));
            kPsi0_2 = dtau * rhs6(E0(i,j) + .5*kE0_1,E2(i,j),Psi0(i,j) + .5*kPsi0_1,Psi1(i,j),Psi2(i,j));
            kE0_3   = dtau * rhs5(E0(i,j) + .5*kE0_2,E2(i,j),Psi0(i,j) + .5*kPsi0_2,Psi1(i,j),Psi2(i,j));
            kPsi0_3 = dtau * rhs6(E0(i,j) + .5*kE0_2,E2(i,j),Psi0(i,j) + .5*kPsi0_2,Psi1(i,j),Psi2(i,j));
            kE0_4   = dtau * rhs5(E0(i,j) + kPsi0_3, E2(i,j),Psi0(i,j) + kPsi0_3,   Psi1(i,j),Psi2(i,j));
            kPsi0_4 = dtau * rhs6(E0(i,j) + kE0_3,   E2(i,j),Psi0(i,j) + kPsi0_3,   Psi1(i,j),Psi2(i,j));
            E0(i+1,j)   = E0(i,j)   + (kE0_1   + 2*kE0_2   + 2*kE0_3   + kE0_4   ) / 6;
            Psi0(i+1,j) = Psi0(i,j) + (kPsi0_1 + 2*kPsi0_2 + 2*kPsi0_3 + kPsi0_4 ) / 6;
        end
        fprintf(1, '\b\b\b\b\b\b\b%06.2f%%', (p + M*(o-1)) * 100.0 /M/N );
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
plot(s, E0(tau == time,:), 'g', s, E1(tau == time,:), 'b', s, E2(tau == time,:), 'r')
grid on
legend('E_0', 'E_1', 'E_2', 'Location', 'southeast')
xlabel('S'); ylabel('E_i');
% print('-f1','E_efficient_RK4','-dpng')

figure(2)
set(gca, 'fontsize', 14)
plot(s, Psi0(tau == time,:), 'g', s, Psi1(tau == time,:), 'b', s, Psi2(tau == time,:), 'r')
grid on
legend('\Psi_0', '\Psi_1', '\Psi_2', 'Location', 'southwest')
xlabel('S'); ylabel('\Psi_i');
% print('-f2','Psi_efficient_RK4','-dpng')

figure(3)
set(gca, 'fontsize', 14)
set(gcf,'renderer','zbuffer');
mesh(s, tau, E2)
set(gca, 'xlim', [0 L], 'ylim', [0 time])
% xlabel('S'); ylabel('\tau'); zlabel('E_2');
% print('-f3','E2_efficient_mesh_RK4','-dpng')

figure(4)
set(gca, 'fontsize', 14)
pcolor(s, tau, E2);
shading interp; 
colorbar
xlabel('S'); ylabel('\tau');
% print('-f4','E2_efficient_map_RK4', '-dpng')
% close all