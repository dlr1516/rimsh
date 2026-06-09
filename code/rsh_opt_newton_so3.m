function [Ropt, data] = rsh_opt_newton_so3(lmax, rsh_src, rsh_dst, x_init, data)
%
%

if nargin < 5
    data = struct();
end

% Parameters of newton optimization in struct data
lmax       = get_opt(data, 'lmax', 5);
maxIter    = get_opt(data, 'maxIter', 100);
tolGrad    = get_opt(data, 'tolGrad', 1e-10);
tolStep    = get_opt(data, 'tolStep', 1e-12);
lambda     = get_opt(data, 'lambda0', 1e-3);
lambdaUp   = get_opt(data, 'lambdaUp', 10);
lambdaDown = get_opt(data, 'lambdaDown', 0.3);
verbose    = get_opt(data, 'verbose', true);

% Checks if data contains the generators of RSH representation of
% Lie algebra so(3) and computes them if they are missing or invalid.
if ~isfield(data, 'Ll') | ~iscell(data.Ll) | length(data.Ll) < lmax+1
    for l = 1:lmax
        L = zeros(2*l+1,2*l+1,3);
        [L(:,:,1), L(:,:,2), L(:,:,3)] = rsh_lie_generator(l);
        data.Ll{l} = L;
    end
end

[r,c] = size(x_init);
if r == 3 & c == 3
  R_init = x_init;
  x_init = log_SO3(R_init);
elseif r * c == 3
  R_init = exp_SO3(x_init);
else
    error('invalid input x_init')
end
  
R = rsh_project_to_so3(R_init); % Initial projection

data.fHist      = zeros(maxIter+1,1);
data.gradHist   = zeros(maxIter+1,1);
data.stepHist   = zeros(maxIter,1);
data.lambdaHist = zeros(maxIter,1);
data.acceptHist = false(maxIter,1);
data.w_cum_Hist = zeros(3, maxIter+1); % Also track the accumulated w

% Get initial metrics
[f, g, H] = rsh_corr_grad_hess(lmax, rsh_src, rsh_dst, x_init);
data.fHist(1) = f;
data.gradHist(1) = norm(g);

current_w = zeros(3,1); % Starting from identity in Lie algebra
% If R0_initial is not identity-based, you might want to initialize current_w
% by logm(R0_initial) and then decomposing it using the generators.
% For simplicity here, we assume R0_initial is the "base" R, and we find the delta relative to it.

if verbose
    fprintf('iter %3d: f = %.16g, ||g|| = %.3e\n', 0, f, norm(g));
end

for k = 1:maxIter
    if norm(g) < tolGrad
        break;
    end

    accepted = false;
    delta = zeros(3,1);

    for trial = 1:20 % Inner damping loop
        A = H - lambda * eye(3);

        if rcond(A) < 1e-14 % Check for near-singularity
            lambda = max(1e-6, lambdaUp * lambda); % Increase damping
            continue;
        end

        delta = -A \ g; % Solve for the step

        % Predicted gain of objective value
        pred = g.' * delta + 0.5 * delta.' * H * delta;

        if pred <= 0 % If predicted gain is not positive, increase damping
            lambda = lambdaUp * lambda;
            continue;
        end

        Rtrial = R * exp_SO3(delta); % Apply the step
        Rtrial = rsh_project_to_so3(Rtrial); % Project to ensure it's a valid R_ell

        [ftrial, ~, ~] = rsh_corr_grad_hess(lmax, rsh_src, rsh_dst, Rtrial); % Evaluate objective at trial point

        ared = ftrial - f; % Actual gain

        rho = ared / pred; % Ratio of actual to predicted gain

        if rho > 1e-4 % If actual gain is positive enough (Levenberg-Marquardt criterion)
            R = Rtrial; % Accept the step
            f = ftrial;
            accepted = true;

            if rho > 0.75 % Very good step, decrease damping
                lambda = max(1e-15, lambdaDown * lambda);
            elseif rho < 0.25 % Not great, increase damping
                lambda = lambdaUp * lambda;
            end

            % If we want to track the total `w`, we need to accumulate these
            % deltas, effectively performing a right-multiplied log-map average
            % For now, we will simply rely on logm at the end for `optimal_w`.
            % Alternatively, one could use Baker-Campbell-Hausdorff if `w` is needed iteratively.
            break;
        else % Step rejected, increase damping
            lambda = lambdaUp * lambda;
        end
    end

    disp('delta | R')
    disp([delta, R])

    data.stepHist(k)   = norm(delta);
    data.lambdaHist(k) = lambda;
    data.acceptHist(k) = accepted;

    if ~accepted % If no step was accepted after many trials, stop
        if verbose
            fprintf('iter %3d: No acceptable step found after %d trials. Stopping.\n', k, trial);
        end
        break;
    end

    [f, g, H] = rsh_corr_grad_hess(lmax, rsh_src, rsh_dst, R); % Recalculate model at new R
    data.fHist(k+1)    = f;
    data.gradHist(k+1) = norm(g);

    if verbose
        fprintf('iter %3d: f = %.16g, ||g|| = %.3e, ||step|| = %.3e, lambda = %.3e\n', ...
            k, f, norm(g), norm(delta), lambda);
    end

    if norm(delta) < tolStep % Check if step size is too small
        break;
    end
end

% Trim history
last = find(data.fHist ~= 0, 1, 'last');
if isempty(last), last = 1; end
data.fHist = data.fHist(1:last);
data.gradHist = data.gradHist(1:last);

last2 = find(data.lambdaHist ~= 0 | data.stepHist ~= 0 | data.acceptHist, 1, 'last');
if isempty(last2), last2 = 0; end
data.stepHist = data.stepHist(1:last2);
data.lambdaHist = data.lambdaHist(1:last2);
data.acceptHist = data.acceptHist(1:last2);

data.R = R;
data.f = f;
data.grad = g;
data.H = H;
Ropt= R;
end