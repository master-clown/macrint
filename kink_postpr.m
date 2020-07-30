function kink_postpr()

%     chl = 1;
%     a = 3*chl;
%     N = 4;
%     mode = 3;
%     
%     if mode == 1
%         load_sfx = 'ubs';
%     else
%         load_sfx = 'ub';
%     end
%     
%     dn_lst = kink_create_out_dirs(false, 1, mode, chl, a, N, load_sfx);
%     
%     kink_lst_fn = strcat(dn_lst(4), "/kink_lst.mat");
%     kink_lst = [];
%     
%     phi_lst = 0:90;
%     for phi = phi_lst
%         kink_mat_fn = strcat(dn_lst(3), sprintf("/phi=%d.mat", phi));
%         if exist(kink_mat_fn, 'file') ~= 2
%             fprintf("%s: '%s' does not exist.\n", datestr(datetime('now')), kink_mat_fn);
%             return;
%         end
% 
%         load(kink_mat_fn, "kink_mat");
%         kink_lst = [ kink_lst, squeeze(kink_mat(4, 1, :)) ];
%     end
%     
%     save(kink_lst_fn, "kink_lst");
    
%------Rotation
%
%---Kinking angle
%     figure;
%     hold on;
%     kink_res_graph(1, 1, 2.5, "red");
%     
%     figure;
%     hold on;
%     kink_res_graph(1, 1, 3, "green");
%     
%     figure;
%     hold on;
%     kink_res_graph(1, 1, 4, "blue");
%     
%     figure;
%     hold on;
%     kink_res_graph(1, 1, 2.5, "red");
%     kink_res_graph(1, 1, 3, "green");
%     kink_res_graph(1, 1, 4, "blue");
%     title("Square cases");
    
    figure;
    hold on;
    plot_kink_sc_rot(1, 'us', "red");
    
    plot_kink(1, 2, 3, "red");
    
    figure;
    hold on;
    plot_kink(1, 3, 3, "blue");
    
    figure;
    hold on;
    plot_kink(1, 2, 3, "red");
    plot_kink(1, 3, 3, "blue");
    title("Rectangular cases");
end

function plot_kink(def_type, mode, a_coef, color)
%
% def_type:
%   1: rotation, 2: y-axis offset.

    chl = 1;
    a = a_coef*chl;
    N = 4;
    
    period_str = sprintf("$%.1fl$", a_coef);
    ls_lst = [ "-", "-.", "--" ];
    legstr_lst = ...
    [   ...
        strcat("Uniaxial ", period_str) ...
        strcat("Shear ", period_str) ...
        strcat("Biaxial ", period_str) ...
    ];
    if(mode ~= 1)
        ls_lst(2) = ls_lst(3);
        legstr_lst(2) = legstr_lst(3);
    end
    
    ylabel("Kinking angle (deg)");
    if def_type == 1
        xlabel("Rotation angle");
        x_lst = 0:90;
        def_type_str = "rotation_test";
    else
        xlabel("Y-axis offset");
        
        x_lst = 0 : 0.1*chl : (a/2);    
        if mode == 3
            x_lst = 2*x_lst;
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
        title("Square lattice, $l_x = 2l_y$", 'Interpreter', 'latex');
        load_sfx = 'ub';
    end

    dn_lst = kink_create_out_dirs(false, def_type, mode, chl, a, N, load_sfx);
    kink_lst_fn = strcat(dn_lst(4), "/kink_lst.mat");
    if exist(kink_lst_fn, 'file') ~= 2
        fprintf("%s: '%s' does not exist.\n", datestr(datetime('now')), kink_lst_fn);
        return;
    end

    load(kink_lst_fn, "kink_lst");

    pl_num = strlength(load_sfx);
    for i_pl = 1:pl_num
        plot(x_lst, kink_lst(i_pl, x_lst+1), ...
            'LineStyle', ls_lst(i_pl), ...
            'Color', color, ...
            'DisplayName', legstr_lst(i_pl));
    end
    
    lg_obj = legend("show");
    set(lg_obj, 'Interpreter', 'latex');
end

function plot_kink_sc_rot(mode, load_str, color)
%
% Plot kink of a single crack being rotated.
%

    chl = 1.0;
    N = 4;
    
    ls_lst = [ "-", "-."];
    legstr_lst = ...
    [   ...
        "Uniaxial SC" ...
        "Shear SC" ...
    ];
    if(mode ~= 1)
        chl = 2.0;
    end
    
    title("Single crack rotation");
    ylabel("Kinking angle (deg)");
    xlabel("Rotation angle");
    x_lst = 0:90;

    kk_sc = zeros(2, length(x_lst));
    
    for(i = 1:length(x_lst))
        kink_mat = eval_kink_single_crack(chl, x_lst(i)*pi/180, 'us');
        kk_sc(:, i) = kink_mat(1, :);
    end
    
    if(contains(load_str, 'u'))
        plot(x_lst, kk_sc(1, x_lst+1), ...
            'LineStyle', ls_lst(1), ...
            'Color', color, ...
            'DisplayName', legstr_lst(1));
    end
    if(contains(load_str, 's'))
        plot(x_lst, kk_sc(2, x_lst+1), ...
            'LineStyle', ls_lst(2), ...
            'Color', color, ...
            'DisplayName', legstr_lst(2));
    end
    
    lg_obj = legend("show");
    set(lg_obj, 'Interpreter', 'latex');
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
    
    res = squeeze(eval_kink(sif_mat));    
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










































