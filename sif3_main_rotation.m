function sif3_main_rotation()

    sif3_main_impl(3.0, 1);
    sif3_main_impl(4.0, 1);
    sif3_main_impl(2.5, 1);
    sif3_main_impl(3.0, 2);
    sif3_main_impl(3.0, 3);
end

function sif3_main_impl(a_coef, mode)

    chl = 1;
    a = a_coef*chl;
    N = 4;
    
    load_sfx = 'y';
    bst_lst = zeros(2, 1);
    bst_lst(2) = 1;

    dn_lst = sif3_create_out_dirs(true, 1, mode, chl, a, N, load_sfx);
    
    phi_lst = 0:90;
    
    for phi = phi_lst
        
        fprintf('%s: calculations for %i degrees...\n', datestr(datetime('now')), phi);
        
        [sif_mat, avg_tr] = kink_process_one(2, mode, phi*pi/180, 0.0, chl, a, N, bst_lst);
        
        sif_mat_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", phi));
        avg_tr_fn = strcat(dn_lst(2), sprintf("/phi=%d.mat", phi));
        
        save(sif_mat_fn, "sif_mat");
        save(avg_tr_fn, "avg_tr");
    end
end












































