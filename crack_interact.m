function [sif_mat, avg_tr] = crack_interact(cc_lst, cp_lst, cl_lst, bst_lst, sif_ind_lst, st_func)
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
%---bst_lst ("boundary stress tensor list"):
%       List of stress tensors of type [ sigma_xx, sigma_xy;  
%                                        sigma_xy, sigma_yy; ] -- boundary 
%       conditions. List goes by the third index.
%
%---sif_ind_lst:
%       List of crack indices for which to evaluate SIFs (heavily affects
%       the performance)
%
%---st_func:
%       A function which returns a list of handles to two ST functions for
%       mode 1 and 2
%
%-Returns:
%---sif_mat:
%       Matrix n-by-4-by-k of SIFs (for modes I and II: opening and sliding) for
%       both ends of every crack, where 'k' equals to number of stress tensors 
%       passed to the function.
%       
%       For i-th crack (i, 1)-th and (i, 2)-th elements are SIF for mode I at 
%       negative and positive ends (-l and +l) of the crack respectively. 
%       The (i, 3)-th and (i, 4)-th elements are for mode II in the similiar way.
%
%---avg_tr:
%       Average tractions on cracks.

%
    
    [st_func_mode1, st_func_mode2] = st_func();
    [sol_lst, rhs_lst] = find_average_tractions(cc_lst, cp_lst, cl_lst, bst_lst, st_func_mode1, st_func_mode2);
    
    fprintf('%s: Starting evaluation of SIFs...\n', datestr(datetime('now')));
    
    N = numel(sif_ind_lst);
    num_rhs = size(rhs_lst, 2);
    sif_neg = zeros(N, 2, num_rhs);
    sif_pos = zeros(N, 2, num_rhs);
    
    % SIFs to get, integrate equations of the system we must
    parfor k = 1:N % parfor
        
        l = cl_lst(sif_ind_lst(k));
        
        for i_rhs = 1:num_rhs
            
            % for negative and positive ends
            f_n = @(s) sqrt((l-s)/(l+s)) * eval_crack_traction(s, sif_ind_lst(k), sol_lst, rhs_lst, ...
                                                               i_rhs, cc_lst, cp_lst, cl_lst, st_func_mode1, st_func_mode2);         
            f_p = @(s) sqrt((l+s)/(l-s)) * eval_crack_traction(s, sif_ind_lst(k), sol_lst, rhs_lst, ...
                                                               i_rhs, cc_lst, cp_lst, cl_lst, st_func_mode1, st_func_mode2);         

            sif_neg(k, :, i_rhs) = integral(f_n, -l, l, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-10) / sqrt(pi*l);
            sif_pos(k, :, i_rhs) = integral(f_p, -l, l, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-10) / sqrt(pi*l);
        end
    end
    
    sif_mat = [ sif_neg(:, 1, :), sif_pos(:, 1, :), sif_neg(:, 2, :), sif_pos(:, 2, :) ];
    avg_tr  = zeros(N, 2, num_rhs);
    for k = 1:num_rhs
        for j = 1:N
            avg_tr(j, 1, k) = sol_lst(2*sif_ind_lst(j)-1, k);
            avg_tr(j, 2, k) = sol_lst(2*sif_ind_lst(j)-0, k);
        end
    end
    
    fprintf('%s: Done.\n', datestr(datetime('now')));
end

function [av_tr_lst, rhs] = find_average_tractions(cc_lst, cp_lst, cl_lst, bst_lst, st_func_mode1, st_func_mode2)
%
%-Desc:
%   Returns the list of average tractions on the cracks. For i-th crack
%   (2*i-1)-th element is normal average traction,
%   (2*i)-th element -- shear one.
%
%   'rhs' is the RHS of the system. In fact this is the external load
%   and is interpreted in similiar way.
%

    N = numel(cl_lst);
    num_rhs = size(bst_lst, 3);
    
    A = eye(2*N, 2*N);
    b = zeros(2*N, num_rhs);
    line_n_od = zeros(1, N);
    line_n_ev = zeros(1, N);
    line_t_od = zeros(1, N);
    line_t_ev = zeros(1, N);
    
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
            b(2*i - 1, i_rhs) = norm' * bst_lst(:,:, i_rhs) * norm;
            b(2*i - 0, i_rhs) = norm' * bst_lst(:,:, i_rhs) * tang;
        end
        
        ind_lst = setdiff(1:N, i);
        
        % fill the row of the matrix
        parfor j = 1:N % parfor
        if j ~= i
            x0_j = cc_lst(j, 1);
            y0_j = cc_lst(j, 2);
            phi_j= cp_lst(j);
            l_j  = cl_lst(j);

            func_nn = @(x,y) norm' * st_func_mode1(x,y, x0_j,y0_j,phi_j,l_j) * norm;
            func_nt = @(x,y) norm' * st_func_mode2(x,y, x0_j,y0_j,phi_j,l_j) * norm;
            func_tn = @(x,y) norm' * st_func_mode1(x,y, x0_j,y0_j,phi_j,l_j) * tang;
            func_tt = @(x,y) norm' * st_func_mode2(x,y, x0_j,y0_j,phi_j,l_j) * tang;
            
            line_n_od(j) = -average_along_crack(func_nn, tcc, tcp, tcl);
            line_n_ev(j) = -average_along_crack(func_nt, tcc, tcp, tcl);
            line_t_od(j) = -average_along_crack(func_tn, tcc, tcp, tcl);
            line_t_ev(j) = -average_along_crack(func_tt, tcc, tcp, tcl);
        end
        end
        
        A(2*i - 1, setdiff(1:2:(2*N), 2*i-1)) = line_n_od(ind_lst);
        A(2*i - 1, setdiff(2:2:(2*N), 2*i-0)) = line_n_ev(ind_lst);
        A(2*i - 0, setdiff(1:2:(2*N), 2*i-1)) = line_t_od(ind_lst);
        A(2*i - 0, setdiff(2:2:(2*N), 2*i-0)) = line_t_ev(ind_lst);
    end
    
    fprintf('%s: Done. Matrix order: %i-by-%i; condition number: %.8f.\n', ...
            datestr(datetime('now')), 2*N, 2*N, cond(A));
    disp('Starting solving the linear system...');
    av_tr_lst = A \ b;
    rhs = b;
    
    fprintf('%s: Done.\n', datestr(datetime('now')));
end

function pt = eval_crack_traction(s, cr_idx, av_tr_lst, rhs_lst, i_rhs, ...
                                  cc_lst, cp_lst, cl_lst, st_func_mode1, st_func_mode2)
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
    pt = [ rhs_lst(2*cr_idx - 1, i_rhs); ...
           rhs_lst(2*cr_idx - 0, i_rhs) ];
    
    tcc = cc_lst(cr_idx,:);                                                              % target crack center (coordinates)
    tcp = cp_lst(cr_idx);
       
    % fill i-th row of the matrix
    for j = 1:N
    if j ~= cr_idx
        x0_j = cc_lst(j, 1);
        y0_j = cc_lst(j, 2);
        phi_j= cp_lst(j);
        l_j  = cl_lst(j);

        av_n = av_tr_lst(2*j - 1, i_rhs);
        av_t = av_tr_lst(2*j - 0, i_rhs);

        term = [ (norm' * st_func_mode1(tcc(1)+cos(tcp)*s, tcc(2)+sin(tcp)*s, x0_j,y0_j,phi_j,l_j) * norm) * av_n + ...
                 (norm' * st_func_mode2(tcc(1)+cos(tcp)*s, tcc(2)+sin(tcp)*s, x0_j,y0_j,phi_j,l_j) * norm) * av_t;  ...
                 (norm' * st_func_mode1(tcc(1)+cos(tcp)*s, tcc(2)+sin(tcp)*s, x0_j,y0_j,phi_j,l_j) * tang) * av_n + ...
                 (norm' * st_func_mode2(tcc(1)+cos(tcp)*s, tcc(2)+sin(tcp)*s, x0_j,y0_j,phi_j,l_j) * tang) * av_t ];
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









































