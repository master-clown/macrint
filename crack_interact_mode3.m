function [sif_mat, avg_tr] = crack_interact_mode3(cc_lst, cp_lst, cl_lst, bc_lst, sif_ind_lst)
%
%-Params:
%	[n here: number of cracks]
%
%---cc_lst ("crack center list"):
%       Matrix n-by-2. i-th row denotes the coordinates of i-th crack's
%       center.
%
%---cp_lst ("crack phi (angle) list"):
%       Vector (n-by-1) of the crack lines' angles (radians).
%
%---cl_lst ("crack length list"):
%       Vector of crack HALF-lengths. 
%!!!    [ATTENTION! HALF-LENGTHS!!!]
%
%---bc_lst ("boundary condition list"):
%       List of column vectors [ sigma_zx;
%                                sigma_zy ]  -- boundary values of antiplane stresses
%
%-Returns:
%---sif_mat:
%       Matrix n-by-2-by-k of SIF mode III values at negative and positive 
%       ends of cracks, where 'k' equals to number of BCs passed to the function.
%
%---avg_tr:
%       Average tractions on cracks.
%

%     test_st_contour();
%     return;
    
    [sol_lst, rhs_lst] = find_average_tractions(cc_lst, cp_lst, cl_lst, bc_lst);
    
    fprintf('%s: Starting evaluation of SIFs...\n', datestr(datetime('now')));
    
    N = numel(sif_ind_lst);
    num_rhs = size(rhs_lst, 2);
    sif_neg = zeros(N, num_rhs);
    sif_pos = zeros(N, num_rhs);
    
    % SIFs to get, integrate equations of the system we must
    for k = 1:N
        
        l = cl_lst(sif_ind_lst(k));
        
        for i_rhs = 1:num_rhs
            
            % for negative and positive ends
            f_n = @(s) sqrt((l-s)/(l+s)) * eval_crack_traction(s, sif_ind_lst(k), sol_lst, rhs_lst, ...
                                                               i_rhs, cc_lst, cp_lst, cl_lst);         
            f_p = @(s) sqrt((l+s)/(l-s)) * eval_crack_traction(s, sif_ind_lst(k), sol_lst, rhs_lst, ...
                                                               i_rhs, cc_lst, cp_lst, cl_lst);         

            sif_neg(k, i_rhs) = integral(f_n, -l, l, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-10) / sqrt(pi*l);
            sif_pos(k, i_rhs) = integral(f_p, -l, l, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-10) / sqrt(pi*l);
        end
    end
    
    sif_mat(:, 1, :) = sif_neg(:, :);
    sif_mat(:, 2, :) = sif_pos(:, :);
    
    avg_tr  = sol_lst;
    
    fprintf('%s: Done.\n', datestr(datetime('now')));
end

function [av_tr_lst, rhs] = find_average_tractions(cc_lst, cp_lst, cl_lst, bc_lst)
%
%-Desc:
%   Returns the list of average tractions on the cracks. 
%
%   'rhs' is the RHS of the system. In fact this is the external load
%   and is interpreted in similiar way.
%

    N = numel(cl_lst);
    num_rhs = size(bc_lst, 2);
    
    A = eye(N, N);
    b = zeros(N, num_rhs);
    line_n = zeros(1, N);
    
    fprintf('%s: Starting matrix assembly...\n', datestr(datetime('now')));
    
    for i = 1:N
    
        tang = [ cos(cp_lst(i));
                 sin(cp_lst(i)) ];
        norm = [-tang(2);
                 tang(1) ];
        tcc = cc_lst(i,:);                                                          % target crack center (coordinates)
        tcp = cp_lst(i);
        tcl = cl_lst(i);
        
        % fill the right hand side
        for i_rhs = 1:num_rhs
            b(i, i_rhs) = dot(norm, bc_lst(:, i_rhs));
        end
        
        ind_lst = setdiff(1:N, i);
        
        % fill the row of the matrix
        for j = 1:N
        if j ~= i
            x0_j = cc_lst(j, 1);
            y0_j = cc_lst(j, 2);
            phi_j= cp_lst(j);
            l_j  = cl_lst(j);

            func = @(x,y) dot(norm, st_antiplane_vec(x,y, x0_j,y0_j,phi_j,l_j));
            
            line_n(j) = -average_along_crack(func, tcc, tcp, tcl);
        end
        end
        
        A(i, ind_lst) = line_n(ind_lst);
    end
    
    fprintf('%s: Done. Matrix order: %i-by-%i; condition number: %.8f.\n', ...
            datestr(datetime('now')), N, N, cond(A));
    disp('Starting solving the linear system...');
    av_tr_lst = A \ b;
    rhs = b;
    
    fprintf('%s: Done.\n', datestr(datetime('now')));
end

function pt = eval_crack_traction(s, cr_idx, av_tr_lst, rhs_lst, i_rhs, ...
                                  cc_lst, cp_lst, cl_lst)
%
%   Evaluate effective traction on the (cr_idx)-th crack at it's 
%   point 's' (natural parameter, which is zero at the crack's 
%   center). 
%

    N = numel(cl_lst);

    tang = [ cos(cp_lst(cr_idx));
             sin(cp_lst(cr_idx)) ];
    norm = [-tang(2);
             tang(1) ];
    pt = rhs_lst(cr_idx, i_rhs);
    
    tcc = cc_lst(cr_idx,:);                                                              % target crack center (coordinates)
    tcp = cp_lst(cr_idx);
       
    % fill i-th row of the matrix
    for j = 1:N
    if j ~= cr_idx
        x0_j = cc_lst(j, 1);
        y0_j = cc_lst(j, 2);
        phi_j= cp_lst(j);
        l_j  = cl_lst(j);

        av_tr = av_tr_lst(j, i_rhs);

        term = dot(norm, st_antiplane_vec(tcc(1)+cos(tcp)*s, tcc(2)+sin(tcp)*s, x0_j,y0_j,phi_j,l_j)) * av_tr;
        pt = pt + term;
    end
    end
    
    return;
end

function val = average_along_crack(func, cc, cp, cl)
%
%-Desc:
%   Averages the func(x, y) along (cr_idx)-th crack, that is, integrates it and
%   divides on the crack length.
%
    cphi = cos(cp);
    sphi = sin(cp);
    
    f = @(s) func(cc(1) + cphi*s, cc(2) + sphi*s);
    
    val = integral(f, -cl, cl, 'ArrayValued', true, ...
                   'RelTol',1e-6,'AbsTol',1e-10) / (2*cl);                          % to pass to the function only scalars, not packs of args
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

function sv = st_antiplane_vec(x_glob, y_glob, x0, y0, phi0, len)
    
    szx = st_antiplane_zx(x_glob, y_glob, x0, y0, phi0, len);
    szy = st_antiplane_zy(x_glob, y_glob, x0, y0, phi0, len);
    
    c = cos(phi0);
    s = sin(phi0);
    
    sv = [ c*szx - s*szy;
           s*szx + c*szy ];
end

function val = st_antiplane_zx(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);    
    minus = sign(x.*y);
    x = abs(x);
    y = abs(y);
    
    M = sqrt( ((x - len).^2 + y.^2).*((x + len).^2 + y.^2) );
    A = x.^2 - y.^2 - len^2;
    
    res = minus.*((y.*sqrt( M + A ) - x.*sqrt( M - A )) ./ (sqrt(2).*M));
    val = real(res);
end

function val = st_antiplane_zy(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    x = abs(x);
    y = abs(y);
    
    M = sqrt( ((x - len).^2 + y.^2).*((x + len).^2 + y.^2) );
    A = x.^2 - y.^2 - len^2;
    
    res = (x.*sqrt( M + A ) + y.*sqrt( M - A )) ./ (sqrt(2).*M) - 1;
    val = real(res);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST FUNCTION
% Run it to test level curves of stress fields (st_ functions).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_st_contour()
    
    x0 = 6.0;
    y0 = 3.0;
    ph = 45*pi/180;
    l = 1.0;
    
    mesh = linspace(-10, 10, 2001);
    mesh_size = size(mesh, 2);
    [X, Y] = meshgrid(mesh, mesh);
    
    ZX = st_antiplane_zx(X, Y, x0, y0, ph, l);
    ZY = st_antiplane_zy(X, Y, x0, y0, ph, l);
    ZX(isnan(ZX)) = 0;
    ZY(isnan(ZY)) = 0;
    
%     Y = -Y;
%     for ix = 1:mesh_size
%         for iy = 1:mesh_size
%             
%             ZX(ix, iy) = st_antiplane_zx(X(ix, iy), Y(ix, iy), x0, y0, ph, l);
%             ZY(ix, iy) = st_antiplane_zy(X(ix, iy), Y(ix, iy), x0, y0, ph, l);
%         end
%     end

    title_lst = { '\sigma_{zx}', '\sigma_{zy}' };
        
    subplot(1, 2, 1);
    contour(X, Y, ZX, [ -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6 ], ...
           'ShowText','on');
    title(title_lst(1));
    xlabel('x');
    ylabel('y');
    axis('equal');
    
    subplot(1, 2, 2);
    contour(X, Y, ZY, [ -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6 ], ...
           'ShowText','on');
    title(title_lst(2));
    xlabel('x');
    ylabel('y');
    axis('equal');
end









































