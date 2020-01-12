%Generacja drugiej harmonicznej o niskiej intesywnosci 
%na podstawie parametrów z 
%Boris P. Antonyuk, Light-Driven Alignment, str. 118–123, Springer 2009.
%
% Autor: Sylwia Majchrowska
% 4 kwietnia 2017r.

clear all;
clc;
format long e

M = 400;
N = 400;
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
E2(:,1) = 0.05;  % E2(0,tau)
Psi2(:,1) = 0;  % Psi2(0,tau)

fprintf(1, '\nStart...       ');
tic
for j = 1:N
    for i = 1:M
        if(j ~= N)
            E2(i,j+1) = E2(i,j) + ...
                dS * (E0(i,j) * sin(Psi0(i,j) - Psi2(i,j)));
            Psi2(i,j+1) = Psi2(i,j) + ...
                dS * (-E0(i,j)/E2(i,j) * cos(Psi0(i,j) - Psi2(i,j)));
        end
        if(i ~= M)
            E0(i+1,j) = E0(i,j) + ...
                dtau * (-E2(i,j)^2 *(E0(i,j) - cos(Psi0(i,j) - ...
                        Psi2(i,j))));
            Psi0(i+1,j) = Psi0(i,j) + ...
                dtau * (-E2(i,j)^2/E0(i,j) * sin(Psi0(i,j) - Psi2(i,j)));
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
% print('-f1','E_week_EM','-dpng')

figure(2)
set(gca, 'fontsize', 14)
plot(s, Psi0(M,:), 'g', s, Psi2(M,:), 'r')
grid on
legend('\Psi_0', '\Psi_2', 'Location', 'northeast')
xlabel('S'); ylabel('\Psi_i');
% print('-f2','Psi_week_EM','-dpng')

figure(3)
set(gca, 'fontsize', 14)
set(gcf,'renderer','zbuffer');
mesh(s, tau, E2)
set(gca, 'xlim', [0 L], 'ylim', [0 time])
xlabel('S'); ylabel('\tau'); zlabel('E_2');
% print('-f3','E2_week_mesh_EM', '-dpng')

figure(4)
set(gca, 'fontsize', 14)
pcolor(s, tau, E2);
shading interp; 
colorbar
xlabel('S'); ylabel('\tau');
% print('-f4','E2_week_map_EM', '-dpng')
% close all