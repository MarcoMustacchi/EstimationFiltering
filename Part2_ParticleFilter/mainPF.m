close all
clearvars
clc

%% Acquire from user
% Acquire prior mean and variance
mu1 = input('Enter the mean of the prior on the first component of x0 (-5 recommended): ');
mu2 = input('Enter the mean of the prior on the second component of x0 (-5 recommended): ');
mu = [mu1; mu2];
sigma2 = input('Enter the variance of the prior on the components of x0 (10 recommended): ');

% Acquire number of measurements N and number of particles Ns
N = input('Enter the number of measures to be used (50 recommended): ');
Ns = input('Enter the number of particles to be used (1e3 recommended): '); 

%% Model definition
% Define the non-linear state space model
f1 = @(x1, x2, k) 0.5*x1 + 25*x1/(1*x1^2 + 1) + x2*cos(29*k);
f2 = @(x1, x2, k) 0.2*x2 + 0.1*x1/(x2^2+1) + 15*sin(10*cos(2*k));
f = @(x, k) [f1(x(1), x(2), k); f2(x(1), x(2), k)];

h = @(x) [sqrt((x(1)+10)^2 + x(2)^2); sqrt(x(1)^2 + (x(2)-10)^2)];

% Define noise variances
sigma2v = 10;
sigma2w = 1;

% Preallocate arrays for performance
x = zeros(2, N);
y = zeros(2, N);

% Generate initial condition x(0) and compute x(1)
x0 = mu + sqrt(sigma2) * randn(2, 1);
x(:, 1) = f(x0, 1) + sqrt(sigma2v)*randn(2, 1);

% Generate the state and output values
for i = 1 : N
    x(:,i+1) = f(x(:,i),i+1) + sqrt(sigma2v)*randn(2, 1);
    y(:,i) = h(x(:,i)) + sqrt(sigma2w)*randn(2, 1);
end
x = x(:, 1:end-1); % Discard last state value in order to have N samples

% Threshold for resampling
NT = Ns/2;

% Generate the first set of particles according to x0 statistics
xk = mu + sqrt(sigma2)*randn(2, Ns);

% Uniform weight at the beginning
wk = ones(1, Ns) ./ Ns;

% Preallocate arrays for performance
xhat = zeros(2, N);
d = zeros(1, N);

for i = 1 : N  
    for j = 1 : Ns
        % q(.|.) = p(x^i_k|x^i_{k-1})
        xk(:,j) = f(xk(:,j), i) + sqrt(sigma2v) * randn(2,1);
        % then weights are proportional only to the likelihood
        wk(j) = wk(j) * mvnpdf(y(:,i)', h(xk(:,j))', sigma2w*eye(2));
    end
    
    % Normalize weights
    t = sum(wk);
    wk = wk ./ t;
    
    % Monte Carlo approximation of the expectation of x under the posterior
    % distribution \pi
    xhat(:,i) = [xk(1,:); xk(2,:)] * wk';
    
    % Plot ''rappresentazione particellare di \pi''
    % figure;
    % scatter3(xk(1,:),xk(2,:),wk);
    
    % Confidence interval calculation
    conf = 0;
    delta = 0.1; % radius increments
    d(i) = 0;
    while conf < 0.95
        conf = 0;
        d(i) = d(i) + delta;
        for j = 1 : Ns
            if norm(xhat(:, i) - xk(:, j)) < d(i)
                conf = conf + wk(j);
            end
        end
    end  
    
    % Resampling
    Neff = 1 / sum(wk.^2);
    if Neff < NT
        c = zeros(Ns, 1);
        for j = 2 : Ns
            c(j) = c(j-1) + wk(j);
        end
        k = 1;
        u = zeros(Ns, 1);
        u(1) = 1/Ns * rand(1);
        for j = 1 : Ns
            u(j) = u(1) + 1/Ns * (j-1);
            while (u(j) > c(k)) && (k < Ns)
                k = k + 1;
            end
            wk(j) = 1/Ns;
        end
    end
end

% Plot results
figure(2);
run('..\Utility\plotProperties.m')
plot(x(1,:), 'LineWidth', 1.5);
plot(xhat(1,:), 'LineWidth', 1.5);
plot(xhat(1,:)+d, 'k-.'); plot(xhat(1,:)-d, 'k-.');
xlabel('Time [k]');
ylabel('Value');
legend('$x_1(k)$','$\hat{x}_1(k)$','Confidence interval at $95\%$',...
    'Interpreter','latex','Location','northwest','FontSize',12);
hold off;
set(gcf, 'Position', get(0, 'Screensize')); % gca not working
exportgraphics(figure(2),'./figure/parameter_estimation_x1.pdf','ContentType','vector')

figure(3);
run('..\Utility\plotProperties.m')
plot(x(2,:), 'LineWidth', 1.5);
plot(xhat(2,:), 'LineWidth', 1.5);
plot(xhat(2,:)+d, 'k-.'); plot(xhat(2,:)-d, 'k-.');
xlabel('Time [k]');
ylabel('Value');
legend('$x_2(k)$','$\hat{x}_2(k)$','Confidence interval at $95\%$',...
    'Interpreter','latex','Location','northwest','FontSize',12);
hold off;
set(gcf, 'Position', get(0, 'Screensize')); % gca not working
exportgraphics(figure(3),'./figure/parameter_estimation_x2.pdf','ContentType','vector')

% var(x(1,10:end)-xhat(1,10:end))
% var(x(2,10:end)-xhat(2,10:end))