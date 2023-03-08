function cod_tables()

    M = [ 1 2 3 4 ];
    chly = 1;
    
    chlx = chly ./ sqrt(M);
    ay = [ 2.4, 2.4, 4.8  ] * chly;
    ax = [ 2, 4, 1 ] .* ay;
    
    load_case = 4;
    vert_cr = 1;
    
    %%%
    bl_col = 1;
    num_cr = vert_cr + 4;
    mat = zeros(length(M)*2, length(ay)*num_cr);
    for i = 1:length(ay)
        
        bl_row = 1;
        for m = 1:length(M)
            res = cod_func(vert_cr, ax(i), ay(i), chlx(m), chly, M(m), load_case);
            mat(bl_row:(bl_row+1), bl_col:(bl_col+num_cr-1)) = res;
            
            bl_row = bl_row + 2;
        end
        
        bl_col = bl_col + num_cr;
    end
end

function res = cod_func(vert_cr, ax, ay, chlx, chly, N_cell_h_cr, load_case)

    load_sfx = 'uhbs';    
    N = 3;

    dn_lst = cod_create_out_dirs(0, vert_cr, ax, ay, chlx, chly, N_cell_h_cr, N, load_sfx);
    load(strcat(dn_lst(1), "/", "avg_tr.mat"), "avg_tr");
    
    if(vert_cr == 1)
        hor_start_idx = 12 + N_cell_h_cr;
        cr_ind_lst = [ 2, hor_start_idx + (0:(N_cell_h_cr-1)) ];
        cod_m_fac = diag([chly, chlx*ones(1, N_cell_h_cr)]);
    else
        hor_start_idx =  2 + N_cell_h_cr;
        cr_ind_lst = [ hor_start_idx + (0:(N_cell_h_cr-1)) ];
        cod_m_fac = diag([chlx*ones(1, N_cell_h_cr)]);
    end
    
    avg_tr = squeeze(avg_tr(cr_ind_lst, :, load_case));
    avg_tr = [ avg_tr(:, 2), avg_tr(:, 1) ];                % 1st col: tang, 2nd: normal
    cod_m = pi*cod_m_fac*avg_tr;
    dE_m = eval_dE(vert_cr, cod_m, ax, ay, chlx, chly, N_cell_h_cr);
    
    num_hor_cr = [ 1 1 2 2 ];       
    if(num_hor_cr(N_cell_h_cr) == 1)
        cod_m(2+vert_cr, :) = NaN(1, 2);
    end
    
    if(vert_cr == 1)
        res = [ cod_m(1,:)', cod_m(2:3, :)', dE_m ];
    else
        res = [ cod_m(1:2, :)', dE_m ];
    end
end

function dE_m = eval_dE(vert_cr, cod_m, ax, ay, lx, ly, N_cell_h_cr)

   % transform from (tang, normal) CS to XY
    cod_xy = cod_m;
    
    if(vert_cr == 1)
        cod_xy(1, :) = [ cod_xy(1, 2), -cod_xy(1, 1) ];
    end
    
    n_hor = [ 0, 1 ];
    n_ver = [ 1, 0 ];
    dE_m = vert_cr*(2*ly*( cod_xy(1,:)'*n_ver + n_ver'*cod_xy(1,:)));
    for(i = (1:N_cell_h_cr) + vert_cr)
        dE_m = dE_m + 2*lx*( cod_xy(i,:)'*n_hor + n_hor'*cod_xy(i,:));
    end
    
    dE_m = dE_m / (2*ax*ay);
end









































