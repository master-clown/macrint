function kink_main()

%     A = zeros(2, 2, 2);
%     A(1, 1, 1) = 1;
%     A(2, 1, 1) = 2;
%     A(1, 2, 1) = 3;
%     A(2, 2, 1) = 4;
%     A(1, 1, 2) = 5;
%     A(2, 1, 2) = 6;
%     A(1, 2, 2) = 7;
%     A(2, 2, 2) = 8;
%     
%     lst = [ squeeze(A(1, 1, :)), squeeze(A(2, 1, :)) ];
%     
% %     
% %     save("results/Am.mat", "A");
% %     load("results/Am.mat", "A");

    kink_main_rect();
end

function kink_main_sq()

    chl = 1;
    a = 2.5*chl;
    %a = 3*chl;
    %a = 4*chl;
    N = 4;
    mode = 1;
    
    load_sfx = 'ubs';
    bst_lst = zeros(2, 2, 3);
    bst_lst(2, 2, 1) = 1;
    bst_lst(1, 2, 2) = 1;
    bst_lst(2, 1, 2) = 1;
    bst_lst(1, 1, 3) = 1;
    bst_lst(2, 2, 3) = 1;

    dn_lst = create_out_dirs(mode, chl, a, N, load_sfx);
    res_fn = strcat(dn_lst(4), "/result.mat");
    
    phi_lst = 55:90;
    if exist(res_fn, 'file') == 2
        load(res_fn, "kink_ang_lst");
    else
        kink_ang_lst = [];
    end
    
    for phi = phi_lst
        
        fprintf('%s: calculations for %i degrees...\n', datestr(datetime('now')), phi);
        
        [sif_mat, avg_tr] = kink_process_one(mode, phi*pi/180, chl, a, N, bst_lst);
        kink_ang = eval_kink(sif_mat / sqrt(pi));
        
        kink_ang_lst = [ kink_ang_lst, squeeze(kink_ang(4, 1, :)) ];
        
        sif_mat_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", phi));
        avg_tr_fn = strcat(dn_lst(2), sprintf("/phi=%d.mat", phi));
        kink_ang_fn = strcat(dn_lst(3), sprintf("/phi=%d.mat", phi));
        
        save(sif_mat_fn, "sif_mat");
        save(avg_tr_fn, "avg_tr");
        save(kink_ang_fn, "kink_ang");
        save(res_fn, "kink_ang_lst");
    end
end

function kink_main_rect()

    chl = 1;
    %a = 2.5*chl;
    a = 3*chl;
    %a = 4*chl;
    N = 4;
    mode = 2;
    
    load_sfx = 'ub';
    bst_lst = zeros(2, 2, 2);
    bst_lst(2, 2, 1) = 1;
    bst_lst(1, 1, 2) = 1;
    bst_lst(2, 2, 2) = 1;

    dn_lst = create_out_dirs(mode, chl, a, N, load_sfx);
    res_fn = strcat(dn_lst(4), "/result.mat");
    
    phi_lst = 57:90;
    if exist(res_fn, 'file') == 2
        load(res_fn, "kink_ang_lst");
    else
        kink_ang_lst = [];
    end
    
    for phi = phi_lst
        
        fprintf('%s: calculations for %i degrees...\n', datestr(datetime('now')), phi);
        
        [sif_mat, avg_tr] = kink_process_one(mode, phi*pi/180, chl, a, N, bst_lst);
        kink_ang = eval_kink(sif_mat / sqrt(pi));
        
        kink_ang_lst = [ kink_ang_lst, squeeze(kink_ang(4, 1, :)) ];
        
        sif_mat_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", phi));
        avg_tr_fn = strcat(dn_lst(2), sprintf("/phi=%d.mat", phi));
        kink_ang_fn = strcat(dn_lst(3), sprintf("/phi=%d.mat", phi));
        
        save(sif_mat_fn, "sif_mat");
        save(avg_tr_fn, "avg_tr");
        save(kink_ang_fn, "kink_ang");
        save(res_fn, "kink_ang_lst");
    end
end

function res = eval_kink(sif_mat)

    sif_siz = size(sif_mat);
    res = zeros(sif_siz(1), 2, sif_siz(3));
    
    for k = 1:sif_siz(3)
        res(:, 1, k) = 2*atan((sif_mat(:, 1, k) - sqrt(...
                               sif_mat(:, 1, k).^2 + 8*sif_mat(:, 3, k).^2) ...
                              ) ./ sif_mat(:, 3, k) / 4);
        res(:, 2, k) = 2*atan((sif_mat(:, 2, k) - sqrt(...
                               sif_mat(:, 2, k).^2 + 8*sif_mat(:, 4, k).^2) ...
                              ) ./ sif_mat(:, 4, k) / 4);
    end
    
    res = 180*res/pi;
end

function dn_lst = create_out_dirs(mode, chl, a, N, load_sfx)
%
%-Params:
%---mode:
%       1: square lattice
%       2: rectangular with a = const and horL = 2 * verL.
%
%---load_sfx:
%       u: uniaxial stretch
%       b: biaxial stretch
%       s: shear

    if mode == 1
        dn_lst =  ...
        [ ...
            sprintf("results/sl_kink/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s", chl, a/chl, N, load_sfx) ...
            sprintf("results/sl_kink/avg_tr/chl=%.1f a=%.1fchl N=%d load=%s", chl, a/chl, N, load_sfx) ...
            sprintf("results/sl_kink/kink_ang/chl=%.1f a=%.1fchl N=%d load=%s", chl, a/chl, N, load_sfx) ...
            sprintf("results/sl_kink/kink_rot_dep/chl=%.1f a=%.1fchl N=%d load=%s", chl, a/chl, N, load_sfx) ...
        ];
    elseif mode == 2
        dn_lst =  ...
        [ ...
            sprintf("results/rl_kink/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s", chl, a/chl, N, load_sfx) ...
            sprintf("results/rl_kink/avg_tr/chl=%.1f a=%.1fchl N=%d load=%s", chl, a/chl, N, load_sfx) ...
            sprintf("results/rl_kink/kink_ang/chl=%.1f a=%.1fchl N=%d load=%s", chl, a/chl, N, load_sfx) ...
            sprintf("results/rl_kink/kink_rot_dep/chl=%.1f a=%.1fchl N=%d load=%s", chl, a/chl, N, load_sfx) ...
        ];
    end
    
    for dn = dn_lst
        mkdir(dn);
    end
end











































