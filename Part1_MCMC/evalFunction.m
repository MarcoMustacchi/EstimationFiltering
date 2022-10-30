function f = evalFunction(x, N)
    f = zeros(N, 1);
    for i = 1 : N
        % f(i) = [1 0] * expm([0 1; -x(1) -x(2)] * i/10) * [x(3); 0];
        % f(i) = exp(x(1)*i/(x(2) + i)) + log(x(3))*i/10;
        % f(i) = x(1)*i / (x(2) + i);
        f(i) = exp(x(2)*i/1000)*x(1) + x(3)/(x(2)^2)*(exp(x(2)*i/1000) -x(2)*i/1000 -1);
    end
end