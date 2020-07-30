function kink_main_offsets()

    kink_main_impl(3.0, 1);
    kink_main_impl(4.0, 1);
    kink_main_impl(2.5, 1);
    kink_main_impl(3.0, 2);
    kink_main_impl(3.0, 3);
end

function kink_main_impl(a_coef, mode)

    chl = 1;
    a = a_coef*chl;
    N = 4;
    
    if mode == 1
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

    dn_lst = kink_create_out_dirs(2, mode, chl, a, N, load_sfx);
    res_fn = strcat(dn_lst(4), "/kink_lst.mat");
    
    off_lst = 0 : 0.1*chl : (a/2);    
    if mode == 3
        off_lst = 2*off_lst;
    end
    
    if exist(res_fn, 'file') == 2
        load(res_fn, "kink_lst");
    else
        kink_lst = [];
    end
    
    off_num = length(off_lst);
    for i = 1:off_num
        
        yoff = off_lst(i);
        
        fprintf('%s: step #%d/%d. Calculations for offset %.3f...\n', datestr(datetime('now')), i, off_num, yoff);
        
        [sif_mat, avg_tr] = kink_process_one(mode, 0.0, yoff, chl, a, N, bst_lst);
        kink_mat = eval_kink(sif_mat);
        
        kink_lst = [ kink_lst, squeeze(kink_mat(4, 1, :)) ];
        
        sif_mat_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", yoff));
        avg_tr_fn = strcat(dn_lst(2), sprintf("/yoff=%.3f.mat", yoff));
        kink_mat_fn = strcat(dn_lst(3), sprintf("/yoff=%.3f.mat", yoff));
        
        save(sif_mat_fn, "sif_mat");
        save(avg_tr_fn, "avg_tr");
        save(kink_mat_fn, "kink_mat");
        save(res_fn, "kink_lst");
    end
end









































































