function kink_main_rotation()

%     kink_main_impl(3.0, 1);
%     kink_main_impl(4.0, 1);
%     kink_main_impl(2.5, 1);
%     kink_main_impl(3.0, 2);
%     kink_main_impl(3.0, 3);

    kink_main_impl(2.5, 4);
end

function kink_main_impl(a_coef, mode)

    chl = 1;
    a = a_coef*chl;
    N = 4;
    
    if mode == 1 || mode == 4
        load_sfx = 'ubs';
        bst_lst = zeros(2, 2, 3);
        bst_lst(2, 2, 1) = 1;
        bst_lst(1, 2, 2) = 1;
        bst_lst(2, 1, 2) = 1;
        bst_lst(1, 1, 3) = 1;
        bst_lst(2, 2, 3) = 1;
    else
        load_sfx = 'ub';
        bst_lst = zeros(2, 2, 2);
        bst_lst(2, 2, 1) = 1;
        bst_lst(1, 1, 2) = 1;
        bst_lst(2, 2, 2) = 1;
    end

    dn_lst = kink_create_out_dirs(1, 1, mode, chl, a, N, load_sfx);
    kink_lst_fn = strcat(dn_lst(4), "/kink_lst.mat");
    
    phi_lst = 0:90;
    
    if exist(kink_lst_fn, 'file') == 2
        load(kink_lst_fn, "kink_lst");
    else
        kink_lst = [];
    end
    
    for phi = phi_lst
        
        fprintf('%s: calculations for %i degrees...\n', datestr(datetime('now')), phi);
        
        [sif_mat, avg_tr] = kink_process_one(1, mode, phi*pi/180, 0.0, chl, a, N, bst_lst);
        kink_mat = eval_kink(sif_mat);
        
        kink_lst = [ kink_lst, squeeze(kink_mat(4, 1, :)) ];
        
        sif_mat_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", phi));
        avg_tr_fn = strcat(dn_lst(2), sprintf("/phi=%d.mat", phi));
        kink_mat_fn = strcat(dn_lst(3), sprintf("/phi=%d.mat", phi));
        
        save(sif_mat_fn, "sif_mat");
        save(avg_tr_fn, "avg_tr");
        save(kink_mat_fn, "kink_mat");
        save(kink_lst_fn, "kink_lst");
    end
end












































