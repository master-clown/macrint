function cod_main()

    M = [ 4 ];
    chly = 1;
    
    chlx = chly ./ sqrt(M);
    ay = [ 2.4, 2.4, 4.8 ] * chly;
    ax = [ 2, 4, 1 ] .* ay;
    
    for i = 1:length(ay)
    for m = 1:length(M)
        cod_func(1, ax(i), ay(i), chlx(m), chly, M(m));
    end
    end
    
%     cod_func(4.8, 4.8, 0.5, 1, 4);
end

function cod_func(vert_cr, ax, ay, chlx, chly, N_cell_h_cr)

    load_sfx = 'uhbs';
    bst_lst = zeros(2, 2, 4);
    bst_lst(2, 2, 1) = 1;
    bst_lst(1, 1, 2) = 1;
    bst_lst(1, 1, 3) = 1;
    bst_lst(2, 2, 3) = 1;
    bst_lst(1, 2, 4) = 1;
    bst_lst(2, 1, 4) = 1;
    
    N = 4;

    fprintf('%s: Calculations for ax = %.1fay, M = %d...\n', datestr(datetime('now')), ax/ay, N_cell_h_cr);

    dn_lst = cod_create_out_dirs(1, vert_cr, ax, ay, chlx, chly, N_cell_h_cr, N, load_sfx);
    sif_fn = strcat(dn_lst(1), "/", "sif_mat.mat");
    avg_fn = strcat(dn_lst(1), "/", "avg_tr.mat");

    [sif_mat, avg_tr] = cod_process(1, vert_cr, ax, ay, chlx, chly, N_cell_h_cr, 1, N, bst_lst);

    save(sif_fn, "sif_mat");
    save(avg_fn, "avg_tr");    
end

function cod_func_simple()

    res_dn = "results/COD/data/single horizontal crack";
    sif_fn = strcat(res_dn, "/", "sif_mat.mat");
    avg_fn = strcat(res_dn, "/", "avg_tr.mat");

    cc_lst = [ 0.0, 0.0 ];
    phi_lst = [ 0 ];
    hl_lst = [ 1 ];
    ind_lst = [ 1 ];
    
    load_sfx = 'uhbs';
    bst_lst = zeros(2, 2, 4);
    bst_lst(2, 2, 1) = 1;
    bst_lst(1, 1, 2) = 1;
    bst_lst(1, 1, 3) = 1;
    bst_lst(2, 2, 3) = 1;
    bst_lst(1, 2, 4) = 1;
    bst_lst(2, 1, 4) = 1;
    
    sif_type = 1;
    
    if(sif_type == 1)
        [sif_mat, avg_tr] = crack_interact(cc_lst, phi_lst, hl_lst, bst_lst, ind_lst);
    else
        [sif_mat, avg_tr] = crack_interact_mode3(cc_lst, phi_lst, hl_lst, bst_lst, ind_lst);
    end

    save(sif_fn, "sif_mat");
    save(avg_fn, "avg_tr");    
end










































