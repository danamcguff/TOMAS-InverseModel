function [Sln, n_iter, err] = ...
    DMPS_TNNLS(Amatrix, Bvector, ApproachPara, MaxIter, Weight, Guess)
% Total Non-Negative Least Square method
% Reference:
% Interior-Point Gradient Method for Large-Scale Totally Nonnegative Least
% Squares Problems (Michael Merritt and Yin Zhang, 2005)
% 
% This method solves problems, which can be described as:
% B = A x Model, where B is measured data (with uncertainties), A is n x m
% matrix, and Model is what we are looking for.
% ApproachParameter is "step" - needs to be smaller than 1, usually 0.6 is
% good. Reasonable range seems to be 0.3 - 0.99
% Limit MaxNumIterations to sensible number...
% Depends on complexity of problem... 
% MaxIter: maximum iteration step;
% Uncert: uncertainties of the vector B, not fully function yet.

% Transcribe from Mandy's Inversion code in Igor by Huajun, 05/18/2016
% Calibrate on 5/20/2016;
% Last update on 5/31/2016;
% Add n_iter and err as output on 4/9/2019 by yhuang


if isempty(Weight)
    A_mat = Amatrix;
    B_vec = Bvector;
    Weight_mat = eye(length(B_vec));
else
    Weight_mat = diag(Weight);
    A_mat = Weight_mat*Amatrix;
    B_vec = Weight_mat*Bvector;
end

A_mat_t = A_mat';

Iter = 1; n_iter = 0;
% [n, m] = size(Amatrix);
% Sln = spline(1:n,Bvector,linspace(1,n,m)');
% Sln = ones(size(Amatrix, 2), 1);
Sln = Guess;
% Initialize the solution

while Iter == 1
    % tic;
    % step (1)
    Qk = A_mat_t * A_mat * Sln - A_mat_t * B_vec;
    Dk = Sln./(A_mat_t * A_mat * Sln);
    Pk = -Dk.*Qk;
    Pk(Qk == 0) = 0;
    % step (2)
    alpha_star = -(Pk' * Qk)/(Pk' * A_mat_t * A_mat * Pk);
    % above is ideal step to make; below is limiting the step so we do not
    % get negative values
    Ak_star = alpha_star * ones(length(Sln), 1);
    alpha_hat = -Sln ./ Pk;
    for ii = 1:length(Pk)
        if Pk(ii) < 0
            Ak_star(ii) = min(ApproachPara * alpha_hat(ii), Ak_star(ii));
        end
    end
    % step (3)
    Sln = Sln + Ak_star.*Pk;
    B_vec_current = Amatrix*Sln;
    div = (Weight_mat*(Bvector-B_vec_current)).^2;
    chisquare = sum(div);
    err = chisquare/length(Bvector);
    n_iter = n_iter + 1;
    Iter = (err > 1) & (n_iter < MaxIter);
    % toc;
end
% fprintf('Iteration: %d; Error: %4.3f; \n', n_iter, err);

if norm(Sln, 1)/length(Sln) <= 1e-32
    Sln = zeros(size(Amatrix, 2), 1);
end
