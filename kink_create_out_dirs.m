function dn_lst = kink_create_out_dirs(create_flag, def_type, mode, chl, a, N, load_sfx)
%
%-Params:
%---def_type:
%       1: rotation, 2: y-axis offset, 3: T-stress.
%
%---mode:
%       1: square lattice
%       2: rectangular with period proportional to l and horL = 2 * verL.
%       3: square with horL = 2 * verL.
%       4: a single row of collinear cracks
%
%---load_sfx:
%       u: uniaxial stretch
%       b: biaxial stretch
%       s: shear

    if(def_type == 1)
        def_type_str = "rotation_test";
    elseif(def_type == 2)
        def_type_str = "offset_test";
    else
        def_type_str = "tstress_test";
    end
    
    if mode == 1
        mode_str = "sl_kink";
    elseif mode == 2
        mode_str = "rl_kink";
    elseif mode == 3
        mode_str = "srl_kink";
    elseif mode == 4
        mode_str = "col_kink";
    end

    dn_lst =  ...
    [ ...
        sprintf("results/%s/%s/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s", def_type_str, mode_str, chl, a/chl, N, load_sfx) ...
        sprintf("results/%s/%s/avg_tr/chl=%.1f a=%.1fchl N=%d load=%s", def_type_str, mode_str, chl, a/chl, N, load_sfx) ...
        sprintf("results/%s/%s/kink_mat/chl=%.1f a=%.1fchl N=%d load=%s", def_type_str, mode_str, chl, a/chl, N, load_sfx) ...
        sprintf("results/%s/%s/kink_lst/chl=%.1f a=%.1fchl N=%d load=%s", def_type_str, mode_str, chl, a/chl, N, load_sfx) ...
    ];
    
    if create_flag == true
        for dn = dn_lst
            mkdir(dn);
        end
    end
end