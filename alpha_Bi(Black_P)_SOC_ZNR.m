clear, clc;
close all

% Parameters
width = 200;
a1 = 3.28373;
a2 = 4.63208;
t0 = 1.0;  %% hopping parameter (eV)
t1 = -1.62;
e = 0.396 * t0;
l = 0.1 * t0; %% SOC strength
B = 0.0;
echarge = 1.6e-19;
hbar = 1.05457e-34;
graphenebond = 1.42e-10;  % Adjust as needed
u = 1;
d = -1;

% Define range of k values
k_values = linspace(0, 2*pi, 5000); % Adjust as needed

% Choose a specific index to print the matrix
i_value = 600;  % Choose an appropriate index based on your k_values

% Preallocate array to store eigenvalues
eigenvaluesU = zeros(length(k_values), 4 * width);
eigenvaluesD = zeros(length(k_values), 4 * width);
%%
% Compute eigenvalues for each k
parfor i = 1:length(k_values)
    k = k_values(i);
    
    % Initialize matrix
    fullMatrixU = zeros(4 * width);
    fullMatrixD = zeros(4 * width);

    % Construct Hamiltonian
    f1 = 2*t0*cos(k.*(1.64186));
    f2 = t1;
    f3 = 2*l*sin(a1.*k); %% SOC function

    % Diagonal block
    for j = 0:width-1
        diagonalBlockU = [e+(u*f3),   f1,    0,    0; 
                          f1', -e-(u*f3),    f2,    0;
                          0,      f2,       u*f3,  f1;
                          0,       0,        f1', -(u*f3)];

        diagonalBlockD = [e+(d*f3), f1,   0,    0; 
                          f1', -e-(d*f3), f2,   0;
                          0,    f2,       d*f3, f1;
                          0,     0,      f1', -(d*f3)];

        fullMatrixU(1 + 4 * j:4 + 4 * j, 1 + 4 * j:4 + 4 * j) = diagonalBlockU;
        fullMatrixD(1 + 4 * j:4 + 4 * j, 1 + 4 * j:4 + 4 * j) = diagonalBlockD;
    end

    % Above diagonal block
    aboveDiagonalBlock = [0, 0, 0, 0;
                          0, 0, 0, 0;
                          0, 0, 0, 0;
                          f2, 0, 0, 0];
    belowDiagonalBlock = [0, 0, 0, f2;
                          0, 0, 0, 0;
                          0, 0, 0, 0;
                          0, 0, 0, 0];
    if width >= 2
        for j = 0:width-2
            fullMatrixU(1 + 4 * j:4 + 4 * j, 5 + 4 * j:8 + 4 * j) = aboveDiagonalBlock;
            fullMatrixU(5 + 4 * j:8 + 4 * j,1 + 4 * j:4 + 4 * j) = belowDiagonalBlock;

            fullMatrixD(1 + 4 * j:4 + 4 * j, 5 + 4 * j:8 + 4 * j) = aboveDiagonalBlock;
            fullMatrixD(5 + 4 * j:8 + 4 * j,1 + 4 * j:4 + 4 * j) = belowDiagonalBlock;

        end
    end

    % Compute eigenvalues
    eigenvaluesU(i, :) = sort(real(eig(fullMatrixU)), 'ascend');
    eigenvaluesD(i, :) = sort(real(eig(fullMatrixD)), 'ascend');


     % Save the matrices for the chosen index
    if i == i_value
        writematrix(fullMatrixU, 'fullMatrixU.txt', 'Delimiter', 'tab');
        writematrix(fullMatrixD, 'fullMatrixD.txt', 'Delimiter', 'tab');
% Print the value of f3
    fprintf('The value of f3 for k = %.10f is %.10f\n', k, f3);
    end
end
%%
% Plot eigenvalues
figure;
plot(k_values, (-eigenvaluesU), 'b');
hold on;
plot(k_values, (-eigenvaluesD), 'r');
xlabel('k');
ylabel('Eigenvalues');
title('Eigenvalues vs. k');
grid on;
%%
% Define range of k values of interest
k_min = 2.2;
k_max = 3.5;
indices_of_interest = find(k_values >= k_min & k_values <= k_max);

% Plot eigenvalues of bands with 10 lowest absolute values closest to zero in the specified range
figure;
plot(k_values(indices_of_interest), real(eigenvaluesU(indices_of_interest, 250:650)), 'r');
hold on;
plot(k_values(indices_of_interest), real(eigenvaluesD(indices_of_interest, 250:650)), 'b');

xlabel('k');
ylabel('Eigenvalues');
grid on;

% Set axis limits
xlim([k_min, k_max]);
ylim([-1.5, 1.5]);
pbaspect([1, 1.5, 1]); % [width, height, depth]
