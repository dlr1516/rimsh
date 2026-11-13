function [F] = rsh_rot_to_z_test_axis_rot(l, u)
%RSH_ROT_TO_Z  Real orthogonal change of basis for RSH rotation blocks.
%
%   F = rsh_rot_to_z(l, u)
%
% Inputs
%   l : band/order degree
%   u : 3x1 or 1x3 Rodrigues rotation vector
%
% Output
%   F : (2*l+1)x(2*l+1) real matrix such that
%           M(u) = F * M_z * F.'
%       where M_z is the real block rotation matrix for rotations about z:
%           blocks R(m*theta), m = l,l-1,...,1 and scalar 1 for m=0.
%
% Notes
%   - Requires:
%       rsh_eval(...)
%       rsh_lm_to_index(...)
%       and generators Lx,Ly,Lz for the band l.
%   - Here it is assumed that a helper function
%         [Lx, Ly, Lz] = rsh_generators(l)
%     exists and returns the (2*l+1)x(2*l+1) real generators in the RSH basis.
%   - The column ordering of F is:
%         [sqrt(2)*a_l, sqrt(2)*b_l, sqrt(2)*a_{l-1}, sqrt(2)*b_{l-1}, ...
%          sqrt(2)*a_1, sqrt(2)*b_1, a_0]
%
%   with v_m = a_m + i b_m and U v_m = i m theta v_m.

    u = u(:);
    if numel(u) ~= 3
        error('u must be a 3-vector.');
    end

    % Band dimension
    n = 2*l + 1;

    % Trivial case l = 0
    if l == 0
        F = 1;
        return;
    end

    % Rotation magnitude
    theta = norm(u);

    % No rotation: F is not unique. Return identity.
    if theta < 1e-14
        F = eye(n);
        return;
    end

    ux = u(1);
    uy = u(2);
    uz = u(3);
    rho = hypot(ux, uy);

    % RSH generators for this band
    [Lx, Ly, Lz] = rsh_lie_generator(l);

    Lu = (ux * Lx + uy * Ly + uz * Lz) / theta;

    [F2, V2, D2] = eig_to_basis(Lu);
    % disp('F2')
    % disp(F2)
    min_val = 100000;

    if rho > 1e-14
        % Original names of variables (L1 > A, L2 > B): 
        %   A = (-uy * Lx + ux * Ly) / rho;
        %   B = (-ux * uz * Lx - uy * uz * Ly + rho^2 * Lz) / (theta * rho);
        L1b = (-uy * Lx + ux * Ly) / rho;
        L2b = (-ux * uz * Lx - uy * uz * Ly + rho^2 * Lz) / (theta * rho);
    else
        % Axis = +/- z
        L1b = Lx;
        if uz >= 0
            L2b = Ly;
        else
            L2b = -Ly;
        end
    end

    for theta = 0:pi/(360*64):2*pi

        % uxt = ux*cos(theta) - uy*sin(theta);
        % uyt = ux*sin(theta) + uy*cos(theta); 
        % 
        % rho = hypot(uxt, uyt);

        Lu = (ux * Lx + uy * Ly + uz * Lz) / theta;

        % Build L1 and L2 defining L_+^(n) = L1 + i L2
        if rho > 1e-14
            % Original names of variables (L1 > A, L2 > B): 
            %   A = (-uy * Lx + ux * Ly) / rho;
            %   B = (-ux * uz * Lx - uy * uz * Ly + rho^2 * Lz) / (theta * rho);
            L1 = cos(theta) * L1b - sin(theta)*L2b;
            L2 = sin(theta) * L1b + cos(theta)*L2b;
        else
            % Axis = +/- z
            L1 = Lx;
            if uz >= 0
                L2 = Ly;
            else
                L2 = -Ly;
            end
        end
    
        % disp('L1 | L2, | Lu')
        % disp([L1, L2, Lu])
        % disp('[L1, L2]')
        % disp(L1 * L2 - L2 * L1)
        % disp('[L2, Lu]')
        % disp(L2 * Lu - Lu * L2)
        % disp('[Lu, L1]')
        % disp(Lu * L1 - L1 * Lu)
    
        %[F2, V2, D2] = eig_to_basis(Lu);
        % disp('F2')
        % disp(F2)
    
        % a0 = v0 from RSH evaluation along axis u
        v_all = rsh_eval(l, u);
        idx0 = sh_lm_to_index(l, -l:l);
        a0 = v_all(idx0);
        a0 = a0(:);
    
        % Normalize a0
        a0 = a0 / norm(a0);
    
        % Storage for a_m, b_m
        a = cell(l+1, 1);
        b = cell(l+1, 1);
    
        a{1} = a0;              % corresponds to m = 0
        b{1} = zeros(n, 1);     % b_0 = 0
    
        % Recurrence for m = 0,1,...,l-1
        for m = 0:(l-1)
            % disp('-------')
            % disp(['m: ' num2str(m)])
            Nm = sqrt((l - m) * (l + m + 1));
    
            a_curr = a{m+1};
            b_curr = b{m+1};
    
            % v_curr = a_curr + i * b_curr;
            % disp('Lu * v_curr | i m * v_curr')
            % disp([Lu * v_curr, i * m * v_curr])
    
            % Name mapping: am > a_curr, bm > b_curr, 
            % a_next = (A * am - B * bm) / Nm;
            % b_next = (A * bm + B * am) / Nm;
            a_next = (L2 * a_curr - L1 * b_curr) / Nm;
            b_next = (L2 * b_curr + L1 * a_curr) / Nm;
    
            % Check for degeneracy
            sa = norm(a_next);
            sb = norm(b_next);
    
            if sa < 1e-14 || sb < 1e-14
                error('Degenerate recurrence at m = %d. Check generator convention.', m+1);
            end
    
            if b_next' * Lu * a_next > 0
                b_next = -b_next;
            end
    
            a{m+2} = a_next;
            b{m+2} = b_next;
        end
    
        % Assemble F in the order:
        % [sqrt(2)a_l, ..., sqrt(2)a_1, a_0, sqrt(2)b_1, ..., sqrt(2)b_l]
        Ftmp = zeros(n, n);
        Ftmp(:,l+1) = a0;
        for m = l:-1:1
            Ftmp(:,l+1-m) = sqrt(2) * a{m+1};
            Ftmp(:,l+1+m) = sqrt(2) * b{m+1};
        end
        
        %exclude center column due to possible sign difference, 
        % these columns are equals in module by definition
        diff = (Ftmp(:,[1:l,l+2:end]) - F2(:,[1:l,l+2:end]));

        cost = sum(abs(diff), 'all');

        if cost < min_val
            min_val = cost;
            Fopt = Ftmp;
        end
    end
    
    F = Fopt;
    disp('F | F2')
    disp([F,F2])
    fprintf('F cost:%f\n', min_val)
    % disp('F | F2')
    % disp([F, F2])
end