function sif_mat = crack_interact(crack_center_lst, crack_phi_lst,  ... 
                                  crack_len_lst, bst)
%
%-Params:
%	[n here: number of cracks]
%
%---crack_center_lst:
%       Matrix n-by-2. i-th row denotes the coordinates of i-th crack's
%       center.
%
%---crack_phi_lst:
%       Vector (n-by-1) of the crack lines' angles (radians).
%
%---crack_len_lst:
%       Vector of crack HALF-lengths. 
%!!!    [ATTENTION! HALF-LENGTHS!!!]
%
%---bst -- "boundary stress tensor":
%       Stress tensor of type [ sigma_xx, sigma_xy;  
%                               sigma_xy, sigma_yy; ] -- boundary conditions.
%
%-Returns:
%   Matrix n-by-4 of SIFs (for modes I and II: opening and sliding) for
%   both ends of every crack. For i-th crack (i, 1)-th and (i, 2)-th 
%   elements are SIF for mode I at negative and positive ends (-l and +l)
%   of the crack respectively. The (i, 3)-th and (i, 4)-th elements are
%   for mode II in the similiar way.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST FUNCTION
% Run it to test level curves of stress fields (st_ functions).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_st_contour()
    
    x0 =-2.0;
    y0 =-4.0;
    ph = pi/8;
    l = 1.0;
    
    mesh = linspace(-10, 10, 1000);
    mesh_size = size(mesh, 2);
    [X, Y] = meshgrid(mesh, mesh);
    
    Z = zeros(mesh_size, mesh_size, 6);
    Z(:, :, 1) = st_normal_xx(X, Y, x0, y0, ph, l);
    Z(:, :, 2) = st_normal_xy(X, Y, x0, y0, ph, l);
    Z(:, :, 3) = st_normal_yy(X, Y, x0, y0, ph, l);
    Z(:, :, 4) = st_shear_xx (X, Y, x0, y0, ph, l);
    Z(:, :, 5) = st_shear_xy (X, Y, x0, y0, ph, l);
    Z(:, :, 6) = st_shear_yy (X, Y, x0, y0, ph, l);
    
    title_lst = { '\sigma^{n}_{xx}', '\sigma^{n}_{xy}', '\sigma^{n}_{yy}', ...
                  '\sigma^{t}_{xx}', '\sigma^{t}_{xy}', '\sigma^{t}_{yy}' };
        
    for i = 1:6
        subplot(2, 3, i);
        contour(X, Y, Z(:,:,i), [ -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6 ], ...
               'ShowText','on');
        title(title_lst(i));
        xlabel('x');
        ylabel('y');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% THE IMPLEMENTATION OF THE METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = eval_along_crack(func, cr_idx, s)
%
%-Desc:
%   Evaluates the func(x, y) at the point 's' of the (cr_idx)-th crack line.
%
%-Params:
%---s:
%       Natural parameter of the crack line. When equal to zero, points to the
%       center of the crack.
%

    x0  = crack_center_lst(cr_idx, 1);
    y0  = crack_center_lst(cr_idx, 2);
    phi = crack_phi_lst(cr_idx);
    cphi= cos(phi);
    sphi= sin(phi);
    
    val = func(x0 + cphi*s, y0 + sphi*s);
end

function val = average_along_crack(func, cr_idx)
%
%-Desc:
%   Averages the func(x, y) along (cr_idx)-th crack, that is, integrates it and
%   divides on the crack length.
%
% Code dublicated from 'eval_along_crack()' for much efficency.
%

    x0  = crack_center_lst(cr_idx, 1);
    y0  = crack_center_lst(cr_idx, 2);
    l   = crack_len_lst(cr_idx);
    phi = crack_phi_lst(cr_idx);
    cphi= cos(phi);
    sphi= sin(phi);
    
    f = @(s) func(x0 + cphi*s, y0 + sphi*s);
    
    val = integral(f, -l, l, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-10) / (2*l);                          % to pass to the function only scalars, not packs of args
end

function [av_tr_lst, rhs] = find_average_tractions()
%
%-Desc:
%   Returns the list of average tractions on the cracks. For i-th crack
%   (2*i-1)-th element is normal average traction,
%   (2*i)-th element -- shear one.
%
%   'rhs' is the RHS of the system. In fact this is the external load
%   and is interpreted in similiar way.
%

    N = numel(crack_len_lst);
    
    A = eye(2*N, 2*N);
    b = zeros(2*N, 1);
    
    for i = 1:N
    
        tang = [ cos(crack_phi_lst(i));
                 sin(crack_phi_lst(i)) ];
        norm = [-tang(2);
                 tang(1) ];
        
        % fill the right hand side
        b(2*i - 1) = norm' * bst * norm;
        b(2*i - 0) = norm' * bst * tang;
        
        % fill the row of the matrix
        for j = 1:N
            if i ~= j
                x0_j = crack_center_lst(j, 1);
                y0_j = crack_center_lst(j, 2);
                l_j  = crack_len_lst(j);
                phi_j= crack_phi_lst(j);
                
                func_nn = @(x,y) norm' *st_normal_tensor(x,y, x0_j,y0_j,phi_j,l_j)*norm;
                func_nt = @(x,y) norm' *st_shear_tensor (x,y, x0_j,y0_j,phi_j,l_j)*norm;
                func_tn = @(x,y) norm' *st_normal_tensor(x,y, x0_j,y0_j,phi_j,l_j)*tang;
                func_tt = @(x,y) norm' *st_shear_tensor (x,y, x0_j,y0_j,phi_j,l_j)*tang;
                
                A(2*i - 1, 2*j - 1) = -average_along_crack(func_nn, i);
                A(2*i - 1, 2*j - 0) = -average_along_crack(func_nt, i);
                A(2*i - 0, 2*j - 1) = -average_along_crack(func_tn, i);
                A(2*i - 0, 2*j - 0) = -average_along_crack(func_tt, i);
            end
        end
    end
    
    fprintf('%s: Matrix order: %i-by-%i; condition number: %.8f\n', ...
            datestr(datetime('now')), 2*N, 2*N, cond(A));
    
    av_tr_lst = A \ b;
    rhs = b;
end

function pt = eval_crack_traction(s, cr_idx, av_tr_lst, rhs)
%
%   Evaluate effective traction on the (cr_idx)-th crack at it's 
%   point 's' (natural parameter, which is zero at the crack's 
%   center). 
%

    N = numel(crack_len_lst);

    tang = [ cos(crack_phi_lst(cr_idx));
             sin(crack_phi_lst(cr_idx)) ];
    norm = [-tang(2);
             tang(1) ];
    pt = [ rhs(2*cr_idx - 1); ...
           rhs(2*cr_idx - 0) ];
    
    % fill i-th row of the matrix
    for j = 1:N
        if cr_idx ~= j
            x0_j = crack_center_lst(j, 1);
            y0_j = crack_center_lst(j, 2);
            l_j  = crack_len_lst(j);
            phi_j= crack_phi_lst(j);

            func_nn = @(x,y) norm' *st_normal_tensor(x,y, x0_j,y0_j,phi_j,l_j)*norm;
            func_nt = @(x,y) norm' *st_shear_tensor (x,y, x0_j,y0_j,phi_j,l_j)*norm;
            func_tn = @(x,y) norm' *st_normal_tensor(x,y, x0_j,y0_j,phi_j,l_j)*tang;
            func_tt = @(x,y) norm' *st_shear_tensor (x,y, x0_j,y0_j,phi_j,l_j)*tang;
            
            av_n = av_tr_lst(2*j - 1);
            av_t = av_tr_lst(2*j - 0);
            
            term = [ eval_along_crack(func_nn, cr_idx, s)*av_n + ...
                     eval_along_crack(func_nt, cr_idx, s)*av_t;  ...
                     eval_along_crack(func_tn, cr_idx, s)*av_n + ...
                     eval_along_crack(func_tt, cr_idx, s)*av_t ];
            pt = pt + term;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% main
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %test_st_contour();

    [sol, rhs] = find_average_tractions();
    
    N = numel(crack_len_lst);
    sif_mat = zeros(N, 4);
    
    % SIFs to get, integrate equations of the system we must
    for k = 1:N
        
        l = crack_len_lst(k);
        
        f_n = @(s) sqrt((l-s)/(l+s)) * eval_crack_traction(s, k, sol, rhs);         % for negative ends
        f_p = @(s) sqrt((l+s)/(l-s)) * eval_crack_traction(s, k, sol, rhs);         %     positive
        
        sif_mat(k, [1,3]) = integral(f_n, -l, l, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-10) / sqrt(pi*l);
        sif_mat(k, [2,4]) = integral(f_p, -l, l, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-10) / sqrt(pi*l);
    end
end


function [x_loc, y_loc] = crack_lc(x, y, x0, y0, phi0)                              % 'crack local coordinates'
%
% Desc: 
%   Converts coordinates (x, y) from global CS to local CS of the crack,
%   center of which is at (x0, y0) and which rotated on phi0 degrees
%   from the global axis X.
%
    cphi0 = cos(phi0); sphi0 = sin(phi0);
    
    x_loc = cphi0*(x-x0) + sphi0*(y-y0);
    y_loc =-sphi0*(x-x0) + cphi0*(y-y0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STRESS TENSOR FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-Prefix: st ('stress tensor')
%
%-Common params:
%---x0, y0:
%       Coordinates of the center of the crack's line.
%
%---phi0:
%       Angle between the crack's line and axis X.
%
%---len:
%       Length of the crack.
%
%-Tests:
%   See the function 'test_st_contour()'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = st_const_a(x, y, len)                                                % 'alpha'
% Auxiliary constant function
    
    val = (x - len).^2 + y.^2;
end

function val = st_const_b(x, y, len)                                                % 'beta'
% Auxiliary constant function
    
    val = 2*(x.^2 + y.^2 - len.^2);
end

function val = st_const_g(x, y, len)                                                % 'gamma'
% Auxiliary constant function
    
    val = (x + len).^2 + y.^2;
end

function val = st_const_d(x, y, len)                                                % 'delta'
% Auxiliary constant function
    
    val = st_const_b(x, y, len) + ...
          2*sqrt(st_const_a(x, y, len).* ...
                 st_const_g(x, y, len));
end

function val = st_const_I2(x, y, len)
% Auxiliary constant function
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(st_const_d(x,y,len));
    
    val = 4*len^2 * 1 ...
                 ./ (ar + gr + dr) ...
                 ./ dr;
end

function val = st_const_I3(x, y, len)
% Auxiliary constant function
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(st_const_d(x,y,len));
    
    val = 2*len^3 * (gr - ar) ...
                 ./ (ar.*gr.*(dr).^3);
end

function val = st_const_I4(x, y, len)
% Auxiliary constant function
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(st_const_d(x,y,len));
    
    val = 2*len^2 * (gr + ar) ...
                 ./ (ar.*gr.*(dr).^3);
end

function val = st_const_I5(x, y, len)
% Auxiliary constant function
    
    d = st_const_d(x,y,len);
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(d);
        
    val =0.5*len^3 * (d.*(gr.^3 - ar.^3) + 3*ar.*gr.*(ar+gr).^2.*(gr-ar)) ...
                  ./ ((ar.*gr).^3.*(dr).^5);
end

function val = st_const_I6(x, y, len)
% Auxiliary constant function
    
    d = st_const_d(x,y,len);
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(d);
        
    val =0.5*len^2 * (d.*(gr.^3 + ar.^3) + 3*ar.*gr.*(ar+gr).^3) ...
                  ./ ((ar.*gr).^3.*(dr).^5);
end

function val = st_normal_xx(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    
    I2 = st_const_I2(x, y, len);
    I4 = st_const_I4(x, y, len);
    I6 = st_const_I6(x, y, len);
    
    val = I2 - 8*y.^2.*(I4 - y.^2.*I6);    
end

function val = st_normal_xy(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    
    I3 = st_const_I3(x, y, len);
    I4 = st_const_I4(x, y, len);
    I5 = st_const_I5(x, y, len);
    I6 = st_const_I6(x, y, len);
    
    val = 2*(-y.*I3 + x.*y.*I4 + 4*y.^3.*(I5 - x.*I6));
end

function val = st_normal_yy(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    
    I2 = st_const_I2(x, y, len);
    I4 = st_const_I4(x, y, len);
    I6 = st_const_I6(x, y, len);
    
    val = I2 + 4*y.^2.*(I4 - 2*y.^2.*I6);    
end

function st = st_normal_tensor(x_glob, y_glob, x0, y0, phi0, len)
    
    sxx = st_normal_xx(x_glob, y_glob, x0, y0, phi0, len);
    sxy = st_normal_xy(x_glob, y_glob, x0, y0, phi0, len);
    syy = st_normal_yy(x_glob, y_glob, x0, y0, phi0, len);
    
    cphi0 = cos(phi0);
    sphi0 = sin(phi0);
    rot_mat = [ cphi0, sphi0;
               -sphi0, cphi0 ];
    
    st = rot_mat' * [ sxx, sxy;
                      sxy, syy ] * rot_mat;
end

function val = st_shear_xx(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    
    I3 = st_const_I3(x, y, len);
    I4 = st_const_I4(x, y, len);
    I5 = st_const_I5(x, y, len);
    I6 = st_const_I6(x, y, len);
    
    val = 2*(3*y.*I3 - 3*x.*y.*I4 - 4*y.^3.*(I5 - x.*I6));
end

function val = st_shear_xy(x_glob, y_glob, x0, y0, phi0, len)
    
    val = st_normal_xx(x_glob, y_glob, x0, y0, phi0, len);
end

function val = st_shear_yy(x_glob, y_glob, x0, y0, phi0, len)
    
    val = st_normal_xy(x_glob, y_glob, x0, y0, phi0, len);
end

function st = st_shear_tensor(x_glob, y_glob, x0, y0, phi0, len)
    
    sxx = st_shear_xx(x_glob, y_glob, x0, y0, phi0, len);
    sxy = st_shear_xy(x_glob, y_glob, x0, y0, phi0, len);
    syy = st_shear_yy(x_glob, y_glob, x0, y0, phi0, len);
    
    cphi0 = cos(phi0);
    sphi0 = sin(phi0);
    rot_mat = [ cphi0, sphi0;
               -sphi0, cphi0 ];
    
    st = rot_mat' * [ sxx, sxy;
                      sxy, syy ] * rot_mat;
end











































