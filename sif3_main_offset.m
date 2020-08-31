function sif3_main_rotation()

    sif3_main_impl(2.5, 1);
%     sif3_main_impl(0:90, 3.0, 1);
%     sif3_main_impl(19:90, 4.0, 1);
%     sif3_main_impl(0:90, 3.0, 2);
%     sif3_main_impl(0:90, 3.0, 3);
end

function sif3_main_impl(a_coef, mode)

    chl = 1;
    a = a_coef*chl;
    N = 4;
    
    load_sfx = 'y';
    bst_lst = zeros(2, 1);
    bst_lst(2) = 1;

    off_lst = 0 : 0.1*chl : (a/2);    
    if mode == 3
        off_lst = 2*off_lst;
    end
    
    dn_lst = sif3_create_out_dirs(true, 2, mode, chl, a, N, load_sfx);
    
    off_num = length(off_lst);
    for i = 1:off_num
                
        yoff = off_lst(i);
        
        fprintf('%s: step #%d/%d. Calculations for offset %.3f...\n', datestr(datetime('now')), i, off_num, yoff);
        
        [sif_mat, avg_tr] = kink_process_one(2, mode, 0.0, yoff, chl, a, N, bst_lst);
        
        sif_mat_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", yoff));
        avg_tr_fn = strcat(dn_lst(2), sprintf("/yoff=%.3f.mat", yoff));
        
        save(sif_mat_fn, "sif_mat");
        save(avg_tr_fn, "avg_tr");
    end
end












































