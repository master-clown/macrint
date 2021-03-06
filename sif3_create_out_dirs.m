function dn_lst = sif3_create_out_dirs(create_flag, def_type, mode, chl, a, N, load_sfx)
%
%-Params:
%---def_type:
%       1: rotation, 2: y-axis offset
%
%---mode:
%       1: square lattice
%       2: rectangular with period proportional to l and horL = 2 * verL.
%       3: rectangular with constant period and horL = 2 * verL.
%
%---load_sfx:
%       x: ZX antiplane
%       y: ZY antiplane

    if(def_type == 1)
        def_type_str = "rotation_test";
    elseif(def_type == 2)
        def_type_str = "offset_test";
    end
    
    if mode == 1
        mode_str = "square";
    elseif mode == 2
        mode_str = "rect";
    else
        mode_str = "square_orthotropic";
    end

    dn_lst =  ...
    [ ...
        sprintf("results/mode3/%s/%s/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s", def_type_str, mode_str, chl, a/chl, N, load_sfx) ...
        sprintf("results/mode3/%s/%s/avg_tr/chl=%.1f a=%.1fchl N=%d load=%s", def_type_str, mode_str, chl, a/chl, N, load_sfx) ...
    ];
    
    if create_flag == true
        for dn = dn_lst
            mkdir(dn);
        end
    end
end