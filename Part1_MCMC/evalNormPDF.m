function p = evalNormPDF(x, mu, Sigma)
    N = length(x);
    SigmaInv = Sigma\eye(N);
    p = 1/sqrt((2*pi)^N * det(Sigma)) * exp(-0.5 * (x - mu)' * SigmaInv * ...
        (x - mu));
end