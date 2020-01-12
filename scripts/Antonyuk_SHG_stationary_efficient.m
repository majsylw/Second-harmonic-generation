% Rozwiazanie rownania dla obwiedni amplitudy wyindukowanego pola
% - eksponencjalne dazenie do stanu stacjonarnego
% dla przypadku wysokiej intensywnosci SH
% na podstawie:
% Boris P. Antonyuk, Light-Driven Alignment, str. 118–128, Springer 2009
%
% Autor: Sylwia Majchrowska
% 7 kwietnia 2017r.

alpha2 = 1;
A2 = 1;
u = 1;
psi2 = 1;
psi1 = 1;

rhs = @(t, psi) alpha2 * A2^2 * (u * exp(1i * (psi2 - 2*psi1)) - psi);
[t,y] = ode45(rhs, [0 10], 0);
Edc = u * exp(1i * (psi2 - 2*psi1)) .* ones(1, length(t));

figure(1)
set(gca, 'fontsize', 16)
plot(t, imag(y), 'r', t, real(y), 'b', 'LineWidth',2)
hold on
plot(t, imag(Edc), '-.r', t, real(Edc), '-.b', 'LineWidth',2)
grid on
legend('\Psi_0', 'A_0', 'Location', 'best')
xlabel('t'); ylabel('A_0, \Psi_0');
print('-f1','Phi_efficient','-dpng')
