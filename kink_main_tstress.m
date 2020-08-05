function kink_main_tstress()

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
    phi = 1.0*pi/180;
    
    load_sfx = 'ub';
    bst_lst = zeros(2, 2, 1);
    bst_lst(2, 2, 1) = 1;

    dn_lst = kink_create_out_dirs(true, 3, mode, chl, a, N, load_sfx);
    kink_lst_fn = strcat(dn_lst(4), "/kink_lst.mat");
    kink_new_lst_fn = strcat(dn_lst(4), "/kink_new_lst.mat");
    
    lam_lst = -1 : 0.1 : 1;
    
    if exist(kink_lst_fn, 'file') == 2
        load(kink_lst_fn, "kink_lst");
        load(kink_new_lst_fn, "kink_new_lst");
    else
        kink_lst = [];
        kink_new_lst = [];
    end
    
    lam_num = length(lam_lst);
    for i = 1:lam_num
        
        bst_lst(1, 1, 1) = lam_lst(i)*bst_lst(2, 2, 1);
        
        fprintf('%s: step #%d/%d. Calculations for lambda = %.1f...\n', datestr(datetime('now')), i, lam_num, lam_lst(i));
        
        [sif_mat, avg_tr] = kink_process_one(mode, phi, 0.0, chl, a, N, bst_lst);
%         load("temp.mat", "sif_mat", "avg_tr");
        
        
        kink_mat = eval_kink(sif_mat);
        kink_new_mat = eval_new_kink(sif_mat, bst_lst(2, 2, 1), lam_lst(i));
        
        kink_lst = [ kink_lst, squeeze(kink_mat(4, 1, :)) ];
        kink_new_lst = [ kink_new_lst, squeeze(kink_new_mat(4, 1)) ];
        
        sif_mat_fn = strcat(dn_lst(1), sprintf("/lam=%.3f.mat", lam_lst(i)));
        avg_tr_fn = strcat(dn_lst(2), sprintf("/lam=%.3f.mat", lam_lst(i)));
        kink_mat_fn = strcat(dn_lst(3), sprintf("/lam=%.3f.mat", lam_lst(i)));
        
        save(sif_mat_fn, "sif_mat");
        save(avg_tr_fn, "avg_tr");
        save(kink_mat_fn, "kink_mat");
        save(kink_lst_fn, "kink_lst");
        save(kink_new_lst_fn, "kink_new_lst");
    end
end

function kink_mat = eval_new_kink(sif_mat, syy, lam)

    r = 0.01;
    th_lst = (-10 : 0.1 : 10) * pi/180;
    
    function kink = eval_kink_sc(k1, k2)
       
        function stt = stress_tt(th)
            
            c = cos(th/2);
            s = sin(th/2);
            
            stt = (k1*c - 3*k2*s) .* c.^2 / sqrt(2*pi*r) + syy*(lam-1)*sin(th).^2;
        end
        
        stt = stress_tt(th_lst);
        [maxval, argmax_idx] = max(stt);
        
        kink = th_lst(argmax_idx);
    end
    
    sif_siz = size(sif_mat);
    kink_mat = zeros(sif_siz(1), 2);
    
    for(i_cr = 1:sif_siz(1))
       
        kink_mat(i_cr, 1) = eval_kink_sc(sif_mat(i_cr, 1, 1), sif_mat(i_cr, 3, 1));
        kink_mat(i_cr, 2) = eval_kink_sc(sif_mat(i_cr, 2, 1), sif_mat(i_cr, 4, 1));
    end
end









































































