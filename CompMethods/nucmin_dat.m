function [L,Sigma,output]  = nucmin_dat(Omega,b,lambda,T,eta_1, n, p)
%% Proximal gradient descent for solving nuclear norm minimization
%% Omega is the non-missing indices of Sigma;
%% b=Sigma(Omega); T is the maximum number of iterations;
%% eta_1 is the initial step size for the low-rank component L.
L = zeros(n,p); c = 0;
output = zeros(0,0);
Sigma = L;
error = 1 / 2 * norm(Sigma(Omega) - b,2)^2 + lambda * sum(svd(L));
for t = 1:T
    grad_L = zeros(n,p);
    grad_L(Omega) = Sigma(Omega) - b;
    if t > 1
        Delta_L = L - L_prev;
        inprod = Delta_L.*(grad_L - grad_L_prev);
        denom1 = sum(inprod(:));
        num1 = (norm(Delta_L,'fro'))^2; 
        eta1_temp = num1 / denom1;
        if ~isnan(eta1_temp)&&(eta1_temp>=1e-4)
            eta_1=eta1_temp;
        end
    end
    L_prev = L; error_prev = error; 
    Sigma_prev = Sigma;
    accept = false;
    while ~accept
        L_temp = L_prev - eta_1 * grad_L;
        L = SVT(L_temp, lambda * eta_1);
        Sigma = L;
        error = 1 / 2 * norm(Sigma(Omega) - b,2)^2 + lambda * sum(svd(L));
        accept = (error <= error_prev)||(eta_1<1e-4);
        eta_1 = eta_1/2;
    end
    grad_L_prev = grad_L;
    grad_norm = norm(Sigma - Sigma_prev, 'fro') / eta_1 / norm(b,2);
    output = cat(1,output,[error, grad_norm, rank(L)]);
    if grad_norm < 1e-8
        break
    end
end
end