close all
clearvars
clc

%% Acquire from user
% Acquire parameters from user
theta1 = input('Enter the value of theta1 (200 recommended): ');
theta2 = input('Enter the value of theta2 (-15 recommended): ');
theta3 = input('Enter the value of theta3 (250 recommended): ');

theta = [theta1; theta2; theta3];
m = length(theta);

% Acquire number of Markov chain iterations from user
M = input('Enter the number of iterations of the Markov chain to be performed (1e5 recommended): ');

%% a priori
thetaMean = zeros(m, 1);
thetaSigma = 1e6 * eye(m);
thetaSigmaInv = thetaSigma \ eye(m);

N = 1000;       % number of measurements
sigma2w = 10;   % noise variance

burnIn = 1e3;   % Number of burn-in samples

%% Generate and Plot outputs
% Generate clean output
f = evalFunction(theta, N);

% Generate noisy output
y = f + sqrt(sigma2w)*randn(N, 1);

% Plot both outputs
figure(1);
run('..\Utility\plotProperties.m')
plot(f, 'LineWidth', 1.5);
plot(y, 'LineWidth', 1.0);
xlabel('Measurement [i]')
ylabel('Value')
legend('Non-noisy outputs $f_i(\theta)$','Measurements $y_i$','Interpreter','latex');
hold off;
set(gcf, 'Position', get(0, 'Screensize')); % gca not working
exportgraphics(figure(1),'./figure/real_output_and_noise_output.pdf','ContentType','vector')

% Preallocate arrays for performance
x = zeros(m, M);
alpha = zeros(M, 1); % Acceptance probability

% Random walk covariance matrix during the burn-in
SigmaRW = diag([100, 10, 100]);
% SigmaRW = 2.3e-3*SigmaRW; % uncomment this line to see burn-in + line 58

%% Run the Markov chain
for i = 1 : M

    % higher variance during burn-in phase
    % comment this line to see burn-in + line 52
    if i == burnIn
        SigmaRW = 2.3e-3*SigmaRW; % Random walk covariance matrix after the burn-in
    end

    % New proposed sample (random walk)
    c = x(:, i) + chol(SigmaRW, 'lower') * randn(m, 1);

    f_c = evalFunction(c, N);
    f_x = evalFunction(x(:, i), N);
    e_c = y - f_c;
    e_x = y - f_x;

    % Compute acceptance probability
    p = exp(-0.5 * (1/sigma2w * (e_c' * e_c - e_x' * e_x) + ...
        c' * thetaSigmaInv * c - x(:, i)' * thetaSigmaInv * x(:, i)));

    alpha(i) = min(1, p);

    % Move the chain
    u = rand(1);
    if (u <= alpha(i))
        x(:, i+1) = c;
    else
        x(:, i+1) = x(:, i);
    end
end

% Show mean acceptance probability after the burn-in
prob_di_accettazione_media = mean(alpha(burnIn:end))

x1 = x(1, :);
x2 = x(2, :);
x3 = x(3, :);

figure(2)
set(gcf, 'Position', get(0, 'Screensize')); % gca not working
scatter3(x1,x2,x3);
exportgraphics(figure(2),'./figure/scatter_parameters_no_burn-in.pdf','ContentType','vector')

% Discard the burn-in samples
x1 = x1(burnIn : end);
x2 = x2(burnIn : end);
x3 = x3(burnIn : end);

%% Plots
figure(3);
scatter3(x1,x2,x3);
set(gcf, 'Position', get(0, 'Screensize')); % gca not working
exportgraphics(figure(3),'./figure/scatter_parameters_convergence.pdf','ContentType','vector')

% Plot histograms
figure(4);
subplot(1,2,1);
run('..\Utility\plotProperties.m')
xlabel('Iteration [i]')
ylabel('Value')
set(gcf, 'Position', get(0, 'Screensize')); % gca not working
plot(x1)
subplot(1,2,2);
enhancedHist(x1);
exportgraphics(figure(4),'./figure/histogram_x1.pdf','ContentType','vector')

figure(5);
subplot(1,2,1);
run('..\Utility\plotProperties.m')
xlabel('Iteration [i]')
ylabel('Value')
set(gcf, 'Position', get(0, 'Screensize')); % gca not working
plot(x2)
subplot(1,2,2);
enhancedHist(x2);
exportgraphics(figure(5),'./figure/histogram_x2.pdf','ContentType','vector')

figure(6);
subplot(1,2,1);
run('..\Utility\plotProperties.m')
xlabel('Iteration [i]')
ylabel('Value')
set(gcf, 'Position', get(0, 'Screensize')); % gca not working
plot(x3)
subplot(1,2,2);
enhancedHist(x3);
exportgraphics(figure(6),'./figure/histogram_x3.pdf','ContentType','vector')

%% Monte Carlo approx. of theta expectation under the posterior distribution
theta1hat = mean(x1);
theta2hat = mean(x2);
theta3hat = mean(x3);

thetahat = [theta1hat; theta2hat; theta3hat];

% Predicted output
yhat = evalFunction(thetahat, N);

% figure;
% plot(y - yhat);
% 
% % Variance of the error
% Varianza_errore_uscita = var((y - yhat))
% 
% figure
% plot(yhat); hold; plot(y);
