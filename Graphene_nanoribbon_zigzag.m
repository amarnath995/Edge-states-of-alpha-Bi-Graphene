clear, clc;
format long
%%
N = 50;
Numberofunit = 1;
Nunit = 2 * N;
NTOT = Numberofunit * Nunit;
a = 1.42;
%%
%%%%% Pos1 is matrix which consists of the positions (in the x-y plane)
%%% of 2N carbon atoms in a unit cell of a zigzag GNR %%%%
Pos1 = zeros(NTOT, 2);
for m = 0:Numberofunit - 1
    for n = 0:1
        for i = 1:N
            % Adjust x-position for zigzag configuration
            Pos1(i + n * N + m * 2 * N, 1) = n * sqrt(3) * a;
            % Adjust y-position for zigzag configuration
            Pos1(i + n * N + m * 2 * N, 2) = i * a + m * sqrt(3) * a / 2;
        end
    end
end
%%
% proj2, (b)
%%%%%%%%%%%%% the elements of Disij indicate distance between
%%%%%%%%%%%%% atom i of Pos1 and atom j of Pos2 %%%%%%%%%%%%
Dis = zeros(NTOT, NTOT);
for i = 1:NTOT
    for j = 1:NTOT
        Dis(i, j) = sqrt((Pos1(i, 1) - Pos1(j, 1))^2 + (Pos1(i, 2) - Pos1(j, 2))^2);
    end
end

%%
% proj2, (c)

%%%%%%%%%%%%% the elements of Hamij indicating the hopping parameter %%%%%%
%%%%%%%%%%%%% between atom i and atom j %%%%%%%%%%%
t = -2.7;
Ham = zeros(NTOT, NTOT);
for i = 1:NTOT
    for j = 1:NTOT
        % Considering only nearest neighbors
        if Dis(i, j) <= 1.01 * a && i ~= j
            Ham(i, j) = t;
        end
    end
end

% Apply periodic boundary conditions
for i = 1:Nunit
    Ham(i, i + Nunit) = t;
    Ham(i + Nunit, i) = t;
end

KX = linspace(-pi/(3*a), pi/(3*a), 100);
d = zeros(NTOT, length(KX)); % Correct size for storing eigenvalues
for i = 1:length(KX)
    H = Ham + Ham' * exp(1j * KX(i) * sqrt(3) * a) + Ham * exp(-1j * KX(i) * sqrt(3) * a);
    d(:, i) = sort(real(eig(H))); % Extract eigenvalues and sort them
end

%%
figure(1)
plot(KX, d, 'b', 'linewidth', 2)
grid on
