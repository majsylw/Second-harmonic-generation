function result = diff_fitting_AWSHG(param)

n2w = 1.473859;

M = 2^10;
N = 2^10;
L = (pi^2*8*(1e-14)*1e6)/(1.064e-4*n2w)*param(1)^2*1e2;
time = 100;

aA2w = 0.000125; % [db/m]

aIR2w = (7.81e11*exp(-48.48/0.532))*1e-3;
aUV2w = ((154.2*0.12)/(46.6*0.12+60)*1e-2*exp(4.63/0.532^-4))*1e-3;

k0 = L;
%a2w = (aA2w + aIR2w + aUV2w)/4.343/k0;
a2w = 2.5512/4.343/k0;

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
E2(:,1) = param(2);  % E2(0,tau)
Psi2(:,1) = 0;  % Psi2(0,tau)

rhs1 = @(E0,E2,Psi0,Psi2)    E0 * sin(Psi0 - Psi2) - a2w/2*E2;
rhs2 = @(E0,E2,Psi0,Psi2) -E0/E2 * cos(Psi0 - Psi2);
rhs3 = @(E0,E2,Psi0,Psi2) -E2^2 *(E0 - cos(Psi0 - Psi2));
rhs4 = @(E0,E2,Psi0,Psi2) -E2^2/E0 * sin(Psi0 - Psi2);

for j = 1:N
    for i = 1:M
        if(j ~= N)
            kE2_1   = dS * rhs1(E0(i,j),E2(i,j),           Psi0(i,j),Psi2(i,j));
            kPsi2_1 = dS * rhs2(E0(i,j),E2(i,j),           Psi0(i,j),Psi2(i,j));
            kE2_2   = dS * rhs1(E0(i,j),E2(i,j) + .5*kE2_1,Psi0(i,j),Psi2(i,j) + .5*kPsi2_1);
            kPsi2_2 = dS * rhs2(E0(i,j),E2(i,j) + .5*kE2_1,Psi0(i,j),Psi2(i,j) + .5*kPsi2_1);
            kE2_3   = dS * rhs1(E0(i,j),E2(i,j) + .5*kE2_2,Psi0(i,j),Psi2(i,j) + .5*kPsi2_2);
            kPsi2_3 = dS * rhs2(E0(i,j),E2(i,j) + .5*kE2_2,Psi0(i,j),Psi2(i,j) + .5*kPsi2_2);
            kE2_4   = dS * rhs1(E0(i,j),E2(i,j) + kE2_3,    Psi0(i,j),Psi2(i,j) + kPsi2_3);
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
    end
end

result = E2(:, N).^2./(ones(M,1)*param(1)).^2 * 100;

end