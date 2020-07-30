function kink_postpr()

    kink_res_graph();
end

function kink_res_graph(def_type, mode, a_coef)
%
% def_type:
%   1: rotation, 2: y-axis offset.

    hold on;
    
    ls_lst = [ "-", "--", "-." ];
    clr_lst = [ "red", "blue", "green" ];
    
    chl = 1;
    a = a_coef*chl;
    N = 4;
    
    if def_type == 1
        xlabel("Rotation angle");
        ylabel("Kinking angle");
        phi_lst = 0:90;
        def_type_str = "sl_k------------ink";
    else
        xlabel("Y-axis offset");
        ylabel("Kinking angle (deg)");
        
        off_lst = 0 : 0.1*chl : (a/2);    
        if mode == 3
            off_lst = 2*off_lst;
        end
        
        def_type_str = "offset_test";
    end
    
    if mode == 1
        title("Square lattice");
        load_sfx = 'ubs';
    elseif mode == 2
        title("Rectangular lattice");
        load_sfx = 'ub';
    else
        title("Square lattice, $l_x = 2l_y$");
        load_sfx = 'ub';
    end

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
    
    legend("Uniax $2.5l$", "Shear $2.5l$", "Biax $2.5l$", ...
           "Uniax $3l$", "Shear $3l$", "Biax $3l$", ...
           "Uniax $4l$", "Shear $4l$", "Biax $4l$",'Interpreter','latex');
end

function rc_kink_rot_graph()

    hold on;
    xlabel("Rotation angle");
    ylabel("Kinking angle");
    title("Kinking angle (degrees). Rectangular lattice, $a = 3l$, $l_x = 2l_y$", 'Interpreter','latex');
    
    a_lst = [ 3 ];
    ls_lst = [ "-" "--" ];
    
    chl = 1;
    N = 4;
    load_sfx = 'ub';
    phi_lst = 0:90;
    mode = 3;
    
    if mode == 1
        mode_dn = "sl_kink";
    elseif mode == 2
        mode_dn = "rl_kink";
    elseif mode == 3
        mode_dn = "srl_kink";
    end
    
    kk = zeros(strlength(load_sfx), length(phi_lst));
    kk_sc = zeros(strlength(load_sfx), length(phi_lst));
    
    for i_a = 1:length(a_lst)
        
        a = a_lst(i_a)*chl;

        kink_dn = sprintf("results/%s/kink_ang/chl=%.1f a=%.1fchl N=%d load=%s/", mode_dn, chl, a/chl, N, load_sfx);
        if exist(kink_dn, 'dir') ~= 7
            fprintf("%s: directory '%s' does not exist.\n", kink_dn, datestr(datetime('now')));
            return;
        end
        
        for i = 1:length(phi_lst)
            load(strcat(kink_dn, sprintf("phi=%d.mat", phi_lst(i))), "kink_ang");            
            kk(:, i) = squeeze(kink_ang(4, 1, :));
            
            kink_sc = eval_kink_single_crack((mode~=1 + 1)*chl, phi_lst(i)*pi/180, load_sfx);
            kk_sc(:, i) = kink_sc(1, :);
        end
        
        for i_pl = 1:2
            plot(phi_lst, kk(i_pl, :), 'LineStyle', ls_lst(i_pl), 'Color', "blue");
            plot(phi_lst, kk_sc(i_pl, :), 'LineStyle', ls_lst(i_pl), 'Color', "black");
        end
    end
    
    legend("Uniaxial", "Uniaxial. Single crack", "Biaxial", "Biax. Single crack", 'Interpreter','latex');
end

function res = eval_kink_single_crack(chl, phi, load_str)

    function [ k1, k2 ] = sif_tens_sc(lambda)
    %
    % SIFs for tension of single crack. Assuming syy = 1.0
    %
    %-Params:
    %---lambda:
    %       A coef in sig_{xx}^{\inf} = lambda * sig_{yy}^{\inf}
       
        syy = 1.0;
        
        k1 = 0.5*syy*sqrt(pi*chl)*(1 + lambda - (1-lambda)*cos(pi - 2*phi));
        k2 = 0.5*syy*sqrt(pi*chl)*(             (1-lambda)*sin(pi - 2*phi));
    end

    function [ k1, k2 ] = sif_shear_sc()
    %
    % SIFs for shearing of single crack. Assuming sxy = 1.0
    %
       
        sxy = 1.0;
        
        k1 =-sxy*sqrt(pi*chl)*cos(pi/2 - 2*phi);
        k2 = sxy*sqrt(pi*chl)*sin(pi/2 - 2*phi);
    end

    sif_mat = zeros(1, 4, strlength(load_str));
    load_idx = 1;
    
    if contains(load_str, 'u')
        [ k1, k2 ] = sif_tens_sc(0.0);
        
        sif_mat(1, 1:2, load_idx) = k1;
        sif_mat(1, 3:4, load_idx) = k2;
        
        load_idx = load_idx + 1;
    end    
    if contains(load_str, 's')
        [ k1, k2 ] = sif_shear_sc();
        
        sif_mat(1, 1:2, load_idx) = k1;
        sif_mat(1, 3:4, load_idx) = k2;
        
        load_idx = load_idx + 1;
    end    
    if contains(load_str, 'b')
        [ k1, k2 ] = sif_tens_sc(1.0);
        
        sif_mat(1, 1:2, load_idx) = k1;
        sif_mat(1, 3:4, load_idx) = k2;
        
        load_idx = load_idx + 1;
    end
    
    res = squeeze(eval_kink(sif_mat));    
end

function sif_to_kink(chl, a, N, load_sfx, phi_lst)

    sif_dn = sprintf("results/sl_kink/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s/", chl, a/chl, N, load_sfx);
    if exist(sif_dn, 'dir') ~= 7
        fprintf("%s: directory '%s' does not exist.\n", sif_dn, datestr(datetime('now')));
        return;
    end

end

function stresskink_rot_graph()

    function ss = eval_stress(chl, sif, kink_ang)
    % sif: 
    %   2-vec with SIF I and II mode
        
        ang = kink_ang*pi/180;
        
        ss = sif(1)*(3*cos(ang/2) + cos(3*ang/2)) / 4 / sqrt(pi*chl) - ...
             3*sif(2)*(sin(ang/2) + sin(3*ang/2)) / 4 / sqrt(pi*chl);
    end
    
    idx_lst = [ 4 2 9 10 ];                                                 % 4 2 9 10 | 2 9 10 1 13
    tip_lst = [ 1 2 2  1 ];                                                 % 1 2 2  1 | 1 1  2 1  2
    chl_lst = ones(length(idx_lst));
    
    case_lst = [ "Uniaxial load. $a = 3l$" ...
                 "Shear load. $a = 3l$" ...
                 "Biaxial load. $a = 3l$" ];
    ls_lst = [ "-k" "--k" ":k" "-.k" ];
    lbl_lst = [ "(1, 0) left"   "(-1, 0) right" ...
                "(0, 1) bottom" "(0,-1) top" ];
    
    chl = 1;
    a = 3*chl;
    N = 4;
    load_sfx = 'ubs';
    phi_lst = 0:90;
    mode = 1;
    
    if mode ~= 1
        chl_lst(1:2) = 2 * chl_lst(1:2);        % all horizontal are twice longer
    end
    
    sif_dn = sprintf("results/sl_kink/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s/", chl, a/chl, N, load_sfx);
    if exist(sif_dn, 'dir') ~= 7
        fprintf("%s: directory '%s' does not exist.\n", sif_dn, datestr(datetime('now')));
        return;
    end
    kink_dn = sprintf("results/sl_kink/kink_ang/chl=%.1f a=%.1fchl N=%d load=%s/", chl, a/chl, N, load_sfx);
    if exist(kink_dn, 'dir') ~= 7
        fprintf("%s: directory '%s' does not exist.\n", kink_dn, datestr(datetime('now')));
        return;
    end
    
    ss = zeros(size(phi_lst));
    kk = zeros(size(phi_lst));
    
    for i_case = 1:3
        load_idx = i_case;       % 1: uniax, 2: shear, 3: biax
        figure();
        
        for i_tip = 1:length(idx_lst)

            for i = 1:length(phi_lst)
                sif_fn = strcat(sif_dn, sprintf("phi=%d.mat", phi_lst(i)));
                kink_fn = strcat(kink_dn, sprintf("phi=%d.mat", phi_lst(i)));

                load(sif_fn, "sif_mat");
                load(kink_fn, "kink_ang");

                kk(i) = kink_ang(idx_lst(i_tip), tip_lst(i_tip), load_idx);
                ss(i) = eval_stress(chl_lst(i_tip), sif_mat(idx_lst(i_tip), [tip_lst(i_tip), tip_lst(i_tip)+2], load_idx), kk(i));
            end

            subplot(2, 1, 1);
            hold on;
            plot(phi_lst, ss, ls_lst(i_tip));


            subplot(2, 1, 2);
            hold on;
            plot(phi_lst, kk, ls_lst(i_tip));
        end
        
        subplot(2, 1, 1);
        legend(lbl_lst, 'Interpreter','latex');
        title(strcat("Relative stress at crack tips. ", case_lst(i_case)), 'Interpreter','latex');

        subplot(2, 1, 2);
        legend(lbl_lst, 'Interpreter','latex');
        title(strcat("Kink at crack tips. ", case_lst(i_case)), 'Interpreter','latex');
    end
end

function ERRk1_rot_graph()

    function err = eval_ERR(sif)
    % sif: 
    %   2-vec with SIF I and II mode
        
        err = sif(1)^2 + sif(2)^2;
    end
    
    idx_lst = [ 4 2 9 10 ];                                                 % 4 2 9 10 | 2 9 10 1 13
    tip_lst = [ 1 2 2  1 ];                                                 % 1 2 2  1 | 1 1  2 1  2
    chl_lst = ones(length(idx_lst));
    
    case_lst = [ "Uniaxial load. $a = 3l$" ...
                 "Shear load. $a = 3l$" ...
                 "Biaxial load. $a = 3l$" ];
    ls_lst = [ "-k" "--k" ":k" "-.k" ];
    lbl_lst = [ "(1, 0) left"   "(-1, 0) right" ...
                "(0, 1) bottom" "(0,-1) top" ];
    
    chl = 1;
    a = 3*chl;
    N = 4;
    load_sfx = 'ubs';
    phi_lst = 0:90;
    mode = 1;
    
    if mode ~= 1
        chl_lst(1:2) = 2 * chl_lst(1:2);        % all horizontal are twice longer
    end

    sif_dn = sprintf("results/sl_kink/sif_mat/chl=%.1f a=%.1fchl N=%d load=%s/", chl, a/chl, N, load_sfx);
    if exist(sif_dn, 'dir') ~= 7
        fprintf("%s: directory '%s' does not exist.\n", sif_dn, datestr(datetime('now')));
        return;
    end
    
    err = zeros(size(phi_lst));
    k1 = zeros(size(phi_lst));
    
    for i_case = 1:3
        load_idx = i_case;       % 1: uniax, 2: shear, 3: biax
        figure();

        for i_tip = 1:length(idx_lst)

            for i = 1:length(phi_lst)
                sif_fn = strcat(sif_dn, sprintf("phi=%d.mat", phi_lst(i)));

                load(sif_fn, "sif_mat");

                err(i) = eval_ERR(sif_mat(idx_lst(i_tip), [tip_lst(i_tip), tip_lst(i_tip)+2], load_idx));
                k1(i) = sif_mat(idx_lst(i_tip), tip_lst(i_tip), load_idx);
            end

            subplot(2, 1, 1);
            hold on;
            plot(phi_lst, err, ls_lst(i_tip));
            
            subplot(2, 1, 2);
            hold on;
            plot(phi_lst, k1, ls_lst(i_tip));
        end

        subplot(2, 1, 1);
        legend(lbl_lst, 'Interpreter','latex');
        title(strcat("ERR. ", case_lst(i_case)), 'Interpreter','latex');

        subplot(2, 1, 2);
        legend(lbl_lst, 'Interpreter','latex');
        title(strcat("SIF I. ", case_lst(i_case)), 'Interpreter','latex');
    end
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

%     function res = eval_kink(sif_mat)
% 
%         sif_siz = size(sif_mat);
%         res = zeros(sif_siz(1), 2, sif_mat(3));
% 
%         num_cr = size(sif_mat, 1);
%         for k = 1:sif_siz(3)
%             for icr = 1:num_cr
%                 ms_bef_f = sif_mat(icr, 3, k) < 0;
%                 ms_bef_s = sif_mat(icr, 1, k) < 0;
%                 res(icr, 1, k) = (-1)^ms_bef_f* ...
%                                  2*atan((sif_mat(icr, 1, k) - (-1)^ms_bef_s*sqrt(...
%                                          sif_mat(icr, 1, k)^2 + 8*sif_mat(icr, 3, k)^2) ...
%                                         ) / sif_mat(icr, 3, k) / 4);
% %                 res(icr, 2, k) = 2*atan((sif_mat(icr, 2, k) - sqrt(...
% %                                        sif_mat(icr, 2, k).^2 + 8*sif_mat(icr, 4, k).^2) ...
% %                                       ) ./ sif_mat(icr, 4, k) / 4);
%             end
%         end
% 
%         res = 180*res/pi;
%     end

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









































