function dn_lst = cod_create_out_dirs(create_flag, vert_cr, ax, ay, chlx, chly, N_cell_h_cr, N, load_sfx)
%
%-Params:
%---vert_cr:
%       1: with vertical cracks, 0: without
%---load_sfx:
%       u: uniaxial stretch
%       b: biaxial stretch
%       s: shear

    if(vert_cr == 1)
        dn_lst =  ...
        [ ...                   %dens variation | edge effect
            sprintf("results/COD/data/new/ax=%.1f ay=%.1f lx=%.1f ly=%.1f M=%d N=%d load=%s", ax, ay, chlx, chly, N_cell_h_cr, N, load_sfx) ...
        ];
    else
        dn_lst =  ...
        [ ...                           %dens variation | edge effect
            sprintf("results/COD/data/no_vert/ax=%.1f ay=%.1f lx=%.1f ly=%.1f M=%d N=%d load=%s", ax, ay, chlx, chly, N_cell_h_cr, N, load_sfx) ...
        ];
    end
    
    if create_flag == true
        for dn = dn_lst
            mkdir(dn);
        end
    end
end