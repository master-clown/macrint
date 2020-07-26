function kink_postpr()

    stresskink_rot_graph();
end

function kink_rot_graph()

    hold on;
    
    a_lst = [ 2.5, 3, 4 ];
    ls_lst = [ "-", "--", "-." ];
    clr_lst = [ "red", "blue", "green" ];
    
    for i_a = 1:length(a_lst)
        chl = 1;
        a = a_lst(i_a)*chl;
        N = 4;
        load_sfx = 'ubs';
        phi_lst = 0:90;

        res_fn = sprintf("results/sl_kink/kink_rot_dep/chl=%.1f a=%.1fchl N=%d load=%s/result.mat", chl, a/chl, N, load_sfx);
        if exist(res_fn, 'file') ~= 2
            fprintf("%s: '%s' does not exist.\n", datestr(datetime('now')));
            return;
        end

        load(res_fn, "kink_ang_lst");
        %kink_ang_lst = 180*kink_ang_lst/pi;

        for i_pl = 1:3
            plot(phi_lst, kink_ang_lst(i_pl, phi_lst+1), 'LineStyle', ls_lst(i_a), 'Color', clr_lst(i_pl));
        end
    end
    
    legend("Uniax $2.5l$", "Shear $2.5l$", "Biax $2.5l$", ...
           "Uniax $3l$", "Shear $3l$", "Biax $3l$", ...
           "Uniax $4l$", "Shear $4l$", "Biax $4l$",'Interpreter','latex');
end

function rc_kink_rot_graph()

    hold on;
    
    a_lst = [ 3 ];
    ls_lst = [ "-" "--" ];
    clr_lst = [ "blue" ];
    
    for i_a = 1:length(a_lst)
        chl = 1;
        a = a_lst(i_a)*chl;
        N = 4;
        load_sfx = 'ub';
        phi_lst = 0:56;

        res_fn = sprintf("results/rl_kink/kink_rot_dep/chl=%.1f a=%.1fchl N=%d load=%s/result.mat", chl, a/chl, N, load_sfx);
        if exist(res_fn, 'file') ~= 2
            fprintf("%s: '%s' does not exist.\n", datestr(datetime('now')));
            return;
        end

        load(res_fn, "kink_ang_lst");
        %kink_ang_lst = 180*kink_ang_lst/pi;

        for i_pl = 1:2
            plot(phi_lst, kink_ang_lst(i_pl, phi_lst+1), 'LineStyle', ls_lst(i_pl), 'Color', clr_lst(i_a));
        end
    end
    
    title("Kinking angle. Rectangular lattice, $a = 3l$, $l_x = 2l_y$", 'Interpreter','latex');
    legend("Uniax $3l$", "Biax $3l$",'Interpreter','latex');
end

function stresskink_rot_graph()

    function ss = eval_stress(chl, sif, kink_ang)
    % sif: 
    %   2-vec with SIF I and II mode
        
        ang = kink_ang*pi/180;
        
        ss = sif(1)*(3*cos(ang/2) + cos(3*ang/2)) / 4 / sqrt(pi*chl) - ...
             3*sif(2)*(sin(ang/2) + sin(3*ang/2)) / 4 / sqrt(pi*chl);
    end

    hold on;
    
    idx_lst = [ 4 8 2 9 ];                                                 % 4 8 2 9 | 8 2 9 1 13
    tip_lst = [ 1 2 2 1 ];                                                 % 1 2 2 1 | 1 1 2 1 2
    lbl_lst = [ "Stress. (1, 0) left","Kink. (1, 0) left",  ...
                "Stress. (1, 0) right", "Kink. (1, 0) right", ...
                "Stress. (0, 1) bottom", "Kink. (0, 1) bottom", ...
                "Stress. (0,-1) top", "Kink. (0,-1) top" ];
    ls_lst = [ "-r", "--r", ":r", "-.r", "-b", "--b", ":b", "-.b" ];
    load_idx = 2;       % 1: uniax, 2: shear, 3: biax
    title_txt = "Relative stress and kink at crack tips during biaxial load and $a = 3l$";

    chl = 1;
    a = 3*chl;
    N = 4;
    load_sfx = 'ubs';
    phi_lst = 0:90;
    
    sif_dn = sprintf("results/sl_kink/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s/", chl, a/chl, N, load_sfx);
    if exist(sif_dn, 'dir') ~= 7
        fprintf("%s: directory '%s' does not exist.\n", sif_dn, datestr(datetime('now')));
        return;
    end
    kink_dn = sprintf("results/sl_kink/kink_ang/chl=%.1f a=%.1fchl N=%d load=%s/", chl, a/chl, N, load_sfx);
    if exist(sif_dn, 'dir') ~= 7
        fprintf("%s: directory '%s' does not exist.\n", kink_dn, datestr(datetime('now')));
        return;
    end
    
    ss = zeros(size(phi_lst));
    kk = zeros(size(phi_lst));
    
    for i_tip = 1:length(idx_lst)
        
        for i = 1:length(phi_lst)
            sif_fn = strcat(sif_dn, sprintf("phi=%d.mat", phi_lst(i)));
            kink_fn = strcat(kink_dn, sprintf("phi=%d.mat", phi_lst(i)));

            load(sif_fn, "sif_mat");
            load(kink_fn, "kink_ang");
            
            kk(i) = kink_ang(i_tip, tip_lst(i_tip), load_idx);
            ss(i) = eval_stress(chl, sif_mat(i_tip, [tip_lst(i_tip), tip_lst(i_tip)+2], load_idx), kk(i));
        end
        
        plot(phi_lst, ss, ls_lst(i_tip));
        plot(phi_lst, kk, ls_lst(i_tip+4));
    end
    
    title(title_txt, 'Interpreter','latex');
    legend(lbl_lst, 'Interpreter','latex');
end

function kink_sif_graph()

    chl = 1;
    a = 3*chl;
    N = 4;
    load_sfx = 'ubs';
    
    sif_dn = sprintf("results/sl_kink/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s/", chl, a/chl, N, load_sfx);
    if exist(sif_dn, 'dir') ~= 7
        fprintf("%s: directory '%s' does not exist.\n", datestr(datetime('now')));
        return;
    end
    
    phi_lst = 0:5;
    sif1 = zeros(3, size(phi_lst, 2));
    sif2 = zeros(3, size(phi_lst, 2));
    
    for i = 1:size(phi_lst, 2)
        
        sif_fn = strcat(sif_dn, sprintf("phi=%d.mat", phi_lst(i)));
        load(sif_fn, "sif_mat");
        
        sif1(:, i) = squeeze(sif_mat(4, 1, :));
        sif2(:, i) = squeeze(sif_mat(4, 3, :));        
    end
    
    hold on;
    lin_lst = [ "-", "--", ":" ];
    sym_lst = [ "x", "o" ];
    
    for i = 1:3
        plot(phi_lst,  sif1(i, :), 'Color', 'black', 'LineStyle', lin_lst(i), 'Marker', sym_lst(1));
        plot(phi_lst,  sif2(i, :), 'Color', 'black', 'LineStyle', lin_lst(i), 'Marker', sym_lst(2));
    end
    
    d1 = line(nan, nan, 'Marker', 'none', 'Color', 'black', 'LineStyle', lin_lst(1));
    d2 = line(nan, nan, 'Marker', 'none', 'Color', 'black', 'LineStyle', lin_lst(2));
    d3 = line(nan, nan, 'Marker', 'none', 'Color', 'black', 'LineStyle', lin_lst(3));
    d4 = line(nan, nan, 'Marker', sym_lst(1), 'Color', 'none', 'MarkerEdgeColor', 'black');
    d5 = line(nan, nan, 'Marker', sym_lst(2), 'Color', 'none', 'MarkerEdgeColor', 'black');
    
    legend([ d1, d2, d3, d4, d5 ], ...
           [ "Uniax", "Shear", "Biax", "SIF mode I", "SIF mode II" ]);
end

function res = eval_kink_poly(sif_mat)

    sif_siz = size(sif_mat);
    res = zeros(sif_siz(1), 2, sif_mat(3));

    num_cr = size(sif_mat, 1);
    for k = 1:sif_siz(3)
        for icr = 1:num_cr
            for itp = 1:2
                k1 = sif_mat(icr, itp, k);
                k2 = sif_mat(icr, itp+2, k);
                
                a = k1*k2;
                b = k1^2 - 3*k2^2;
                
                rts = roots([ ...
                              4*a^2 + b^2 ...
                              12*a^2 + 4*b*k2^2 - b^2 ...
                              9*a^2 - 4*b*k2^2 + 4*k2^4 ...
                              -4*k2^4 ...
                            ]);
                rts = 2*acos(sqrt(rts));          % neglected signs after sqrt!
            end
        end
    end

    res = 180*res/pi;
end

function kink_eval_test()

    function res = eval_kink(sif_mat)

        sif_siz = size(sif_mat);
        res = zeros(sif_siz(1), 2, sif_mat(3));

        num_cr = size(sif_mat, 1);
        for k = 1:sif_siz(3)
            for icr = 1:num_cr
                ms_bef_f = sif_mat(icr, 3, k) < 0;
                ms_bef_s = sif_mat(icr, 1, k) < 0;
                res(icr, 1, k) = (-1)^ms_bef_f* ...
                                 2*atan((sif_mat(icr, 1, k) - (-1)^ms_bef_s*sqrt(...
                                         sif_mat(icr, 1, k)^2 + 8*sif_mat(icr, 3, k)^2) ...
                                        ) / sif_mat(icr, 3, k) / 4);
%                 res(icr, 2, k) = 2*atan((sif_mat(icr, 2, k) - sqrt(...
%                                        sif_mat(icr, 2, k).^2 + 8*sif_mat(icr, 4, k).^2) ...
%                                       ) ./ sif_mat(icr, 4, k) / 4);
            end
        end

        res = 180*res/pi;
    end

    chl = 1;
    a = 3*chl;
    N = 4;
    load_sfx = 'ubs';
    
    sif_dn = sprintf("results/sl_kink/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s/", chl, a/chl, N, load_sfx);
    if exist(sif_dn, 'dir') ~= 7
        fprintf("%s: directory '%s' does not exist.\n", datestr(datetime('now')));
        return;
    end
    
    phi_lst = 0:5;
    kink_ang_lst = zeros(3, size(phi_lst, 2));
    
    for i = 1:size(phi_lst, 2)
        
        sif_fn = strcat(sif_dn, sprintf("phi=%d.mat", phi_lst(i)));
        load(sif_fn, "sif_mat");
        
        kink_ang = eval_kink_poly(sif_mat / sqrt(pi));
        
        kink_ang_lst(:, i) = squeeze(kink_ang(4, 1, :));    
    end
    
    plot(phi_lst, kink_ang_lst);
    legend("Uniax", "Shear", "Biax");
end









































