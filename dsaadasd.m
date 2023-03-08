

M = [ 1 ];
chly = 1;

chlx = chly ./ sqrt(M);
ay = [ 2.4 ] * chly;%[ 2.4, 2.4, 4.8  ] * chly;
ax = [ 1 ] .* ay;%[ 2, 4, 1 ] .* ay;


for vert_cr = [1 0]
for iso_nu  = [0.3 0.25]
    mat = zeros(length(M)*2, length(ay)*2, 3);
    
    for load_idx = 1:4
        
        bl_col = 1;
        for i = 1:length(ay)
            
            bl_row = 1;
            for m = 1:length(M)
                res = cod_dE(vert_cr, ax(i), ay(i), chlx(m), chly, M(m), load_idx);
                mat(bl_row:(bl_row+1), bl_col:(bl_col+1), load_idx) = res;
                
                bl_row = bl_row + 2;
            end
            
            bl_col = bl_col + 2;
        end
    end


    % Table format: 5 + 5*8 + 5*8 = 85 symbols width

    if(vert_cr == 1)
        fprintf(sprintf("Eff. elastic constans (rel). Cross-crack confs. Nu = %.2f\n", iso_nu));
    else
        fprintf(sprintf("Eff. elastic constans (rel). Hor.-crack confs. Nu = %.2f\n", iso_nu));
    end
    
%     fprintf("|%3s|%s|%s|%s|\n", "M", strjust(sprintf("%39s", "4x1 Rect"),"center"), ...
%                                      strjust(sprintf("%39s", "2x2 Square"),"center"), ...
%                                      strjust(sprintf("%39s", "2x1 Rect"),"center"));
    fprintf("|  M|");
    for i = 1:1
        fprintf("%8s|%10s|%10s|%10s|%10s|", "Ex", "Ey", "vxy", "vyx", "Gxy");
    end
%     fprintf("|  M|");
%     for i = 1:1
%         fprintf("%8s|%10s|%10s|%10s|%10s|%10s|", "A: LHS", "A: RHS", "B: LHS", "B: RHS", "C: LHS", "C: RHS");
%     end
    fprintf("\n");
    
    mean_perc = zeros(1, 5);
    
    for m = M
        
        dE_2r = mat((2*m-1):(2*m), 1:2, :);
%         dE_4r = mat((2*m-1):(2*m), 3:4, :);
%         dE_2s = mat((2*m-1):(2*m), 5:6, :);
        
%         [Ex1, Ey1, vxy1, vyx1, Gxy1] = eval_rel_eff_consts(dE_4r, iso_nu);
%         [Ex2, Ey2, vxy2, vyx2, Gxy2] = eval_rel_eff_consts(dE_2s, iso_nu);
        [Ex3, Ey3, vxy3, vyx3, Gxy3] = eval_rel_eff_consts(dE_2r, iso_nu);
        
%         l1 = vxy1*Ey1;
%         r1 = vyx1*Ex1;
%         l2 = vxy2*Ey2;
%         r2 = vyx2*Ex2;
        l3 = vxy3*Ey3;
        r3 = vyx3*Ex3;
        
%         perc = abs([ 1-Ex1/Ex2, 1-Ey1/Ey2, 1-vxy1/vxy2, 1-vyx1/vyx2, 1-Gxy1/Gxy2 ]);
%         mean_perc = mean_perc + perc;
        
%         fprintf("    %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n",   Ex1, Ey1, vxy1, vyx1, Gxy1);
%         fprintf("%4i%8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n", m, Ex2, Ey2, vxy2, vyx2, Gxy2);
        %fprintf("   %%%8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n", 1-Ex1/Ex2, 1-Ey1/Ey2, 1-vxy1/vxy2, 1-vyx1/vyx2, 1-Gxy1/Gxy2);
        fprintf("    %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n",   Ex3, Ey3, vxy3, vyx3, Gxy3);     
        %fprintf("   %%%8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n", perc(1), perc(2), perc(3), perc(4), perc(5));  
        %fprintf("    %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n", l1, r1, l2, r2, l3, r3);
    end
    
%     mean_perc = mean_perc / 4;
%     mean = sum(mean_perc) / 4;
    
    %fprintf("\n\n");    
    %fprintf("   %%%8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\  & %8.3f\n", mean_perc(1), mean_perc(2), mean_perc(3), mean_perc(4), mean_perc(5), mean);    
    
    
    rho_x = 4*0.0434028;
    rho_y = rho_x;
    ExNoInt = 1/(1 + 2*pi*rho_y);
    EyNoInt = 1/(1 + 2*pi*rho_x);
    vxyNoInt = ExNoInt;
    vyxNoInt = EyNoInt;
    GxyNoInt = 1/(1 + pi*(rho_x + rho_y)/(1 + iso_nu));
    fprintf("NI: %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n", (1-vert_cr) + vert_cr*ExNoInt, EyNoInt, (1-vert_cr) + vert_cr*vxyNoInt, vyxNoInt, GxyNoInt);
    fprintf("\n\n");
    
end
end

function [Ex, Ey, vxy, vyx, Gxy] = eval_rel_eff_consts(dE_m, iso_nu)
% Those values are RELATIVE!!!
    
    eval_symm = 0;
    method_biax = 0;

    Ex = 1 / (1 + dE_m(1, 1, 2));
    Ey = 1 / (1 + dE_m(2, 2, 1));
    Gxy = 1 / (1 + dE_m(1, 2, 4)/(1+iso_nu)/2);
    
    if(method_biax == 0)
        vxy = (1 - dE_m(2, 2, 2)/iso_nu) * Ex;
        vyx = (1 - dE_m(1, 1, 1)/iso_nu) * Ey;    
    else
        vxy = (iso_nu - dE_m(2, 2, 3) + 1/(Ey) - 1) * Ex / iso_nu;
        vyx = (iso_nu - dE_m(1, 1, 3) + 1/(Ex) - 1) * Ey / iso_nu;
    end
    
    if(eval_symm == 1)
        vxy_old = vxy;
        vyx_old = vyx;
        vxy = 0.5*(vxy_old + vyx_old/(Ey/Ex));
        vyx = 0.5*(vyx_old + vxy_old/(Ex/Ey));
    end
end

function res = cod_dE(vert_cr, ax, ay, chlx, chly, N_cell_h_cr, load_case)

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
    
    res = dE_m;
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
