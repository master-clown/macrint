function kink_postpr()
    
%------Rotation
%
%---Kinking angle
%     figure;
%     hold on;
%     plot_kink_sc_rot(1, 'us', "magenta");
%     plot_kink(1, 1, 2.5, "red");
%     
%     figure;
%     hold on;
%     plot_kink_sc_rot(1, 'us', "magenta");
%     plot_kink(1, 1, 3, "green");
%     
%     figure;
%     hold on;
%     plot_kink_sc_rot(1, 'us', "magenta");
%     plot_kink(1, 1, 4, "blue");
%     
%     figure;
%     hold on;
%     plot_kink(1, 1, 2.5, "red");
%     plot_kink(1, 1, 3, "green");
%     plot_kink(1, 1, 4, "blue");
%     plot_kink_sc_rot(1, 'ubs', "magenta");
%     title("Square cases");
%     
%     figure;
%     hold on;    
%     plot_kink(1, 2, 3, "red");
%     
%     figure;
%     hold on;
%     plot_kink(1, 3, 3, "blue");
%     
%     figure;
%     hold on;
%     plot_kink(1, 2, 3, "red");
%     plot_kink(1, 3, 3, "blue");
%     plot_kink(1, 1, 3, "green", "ub");
%     title("Rectangular cases");
% % 
% % %---Stresses in neighbours
%     plot_stress(1, 1, 2.5, "-");    
%     plot_stress(1, 1, 3, "-");    
%     plot_stress(1, 1, 4, "-");    
%     plot_stress(1, 2, 3, "-");    
%     plot_stress(1, 3, 3, "-");
% %
%---ERR in neighbours
%     plot_err(2, 1, 1, 2.5, "-");    
%     plot_err(1, 1, 1, 3, "-");    
%     plot_err(1, 1, 1, 4, "-");    
%     plot_err(1, 1, 2, 3, "-");    
%     plot_err(1, 1, 3, 3, "-");

%---SIF mode 3 in neighbours
%     figure;
%     plot_sif3(1, 1, 2.5, "-");    
%     
%     figure;
%     plot_sif3(1, 1, 3, "-");    
%     
%     figure;
%     plot_sif3(1, 1, 4, "-");    
%     
%     figure;
%     plot_sif3(1, 2, 3, "-");    
%     
%     figure;
%     plot_sif3(1, 3, 3, "-");

% 
% 
% % %------Y-axis offset
% % %
% % %---Kinking angle
% %     figure;
% %     hold on;
% %     plot_kink(2, 1, 2.5, "red");
% %     
% %     figure;
% %     hold on;
% %     plot_kink(2, 1, 3, "green");
% %     
% %     figure;
% %     hold on;
% %     plot_kink(2, 1, 4, "blue");
% %     
%     figure;
%     hold on;
%     plot_kink_sc_rot(1, 'ubs', "magenta", 0 : 0.1 : 2.0);
%     plot_kink(2, 1, 2.5, "red");
%     plot_kink(2, 1, 3, "green");
%     plot_kink(2, 1, 4, "blue");
%     title("Square cases");
% %     
% %     figure;
% %     hold on;    
% %     plot_kink(2, 2, 3, "red");
% %     
% %     figure;
% %     hold on;
% %     plot_kink(2, 3, 3, "blue");
% %     
%     figure;
%     hold on;
%     plot_kink(2, 1, 3, "green", "ub");
%     plot_kink(2, 2, 3, "red");
%     plot_kink(2, 3, 3, "blue");
%     title("Rectangular cases");
% % 
% %---Stresses in neighbours
%     plot_stress(2, 1, 2.5, "-");    
%     plot_stress(2, 1, 3, "-");    
%     plot_stress(2, 1, 4, "-");    
%     plot_stress(2, 2, 3, "-");    
%     plot_stress(2, 3, 3, "-");
% 
% %---ERR in neighbours
%     plot_err(2, 2, 1, 2.5, "-");    
%     plot_err(2, 2, 1, 3, "-");    
%     plot_err(2, 2, 1, 4, "-");    
%     plot_err(2, 2, 2, 3, "-");    
%     plot_err(2, 2, 3, 3, "-");

% %------T-stress
% %
%     figure;
%     hold on;
%     plot_tstress(1, 2.5, "red");    
%     
%     figure;
%     hold on;
%     plot_tstress(1, 3, "green");    
%     
%     figure;
%     hold on;
%     plot_tstress(1, 4, "blue");    
%     
%     figure;
%     hold on;
%     plot_tstress(2, 3, "black");    
%     
%     figure;
%     hold on;
%     plot_tstress(3, 3, "magenta");
%     %title("T-stress test");

% %---SIF mode 3 in neighbours
%     figure;
%     plot_sif3(2, 1, 2.5, "-");    
    
%     figure;
%     plot_sif3(2, 1, 3, "-");    

%     figure;
%     plot_sif3(2, 1, 4, "-");    

%     figure;
%     plot_sif3(2, 2, 3, "-");    

%     figure;
%     plot_sif3(2, 3, 3, "-");

end

function plot_kink(def_type, mode, a_coef, color, load_str)
%
% def_type:
%   1: rotation, 2: y-axis offset.

    chl = 1;
    a = a_coef*chl;
    N = 4;
    
    ls_lst = [ "-", "-.", "--" ];
    period_str = sprintf("$a = %.1fl$", a_coef);
    if(mode == 3)
        period_str = sprintf("a = %.1fl_x", a_coef);
    end
    legstr_lst = ...
    [   ...
        strcat("Uniaxial ", period_str) ...
        strcat("Shear ", period_str) ...
        strcat("Biaxial ", period_str) ...
    ];

    ylabel("Kinking angle \theta_* [deg]");
    yticks(-180:20:180);
    if(mode ~= 1)
        yticks(-180:10:180);
    end
    
    if def_type == 1
        xlabel("Rotation angle [deg]");
        x_lst = 0:90;
    else
        yticks(-180:10:180);
        xlabel("Y-axis offset [m]");
        
        x_lst = 0 : 0.1*chl : (a/2);    
        if mode == 3
            x_lst = 2*x_lst;
        end
        if mode ~= 1
            yticks(-180:2:180);
        end
        
        xticks(x_lst);
    end
    
    if mode == 1
        title("Square lattice");
        load_sfx = 'ubs';
    elseif mode == 2
        title("Rectangular lattice");
        load_sfx = 'ub';
    else
        title("Square lattice, l_x = 2l_y", 'Interpreter', 'latex');
        load_sfx = 'ub';
    end
    
    if ~exist('load_str','var')
        load_str = load_sfx;
    end
 
    dn_lst = kink_create_out_dirs(false, def_type, mode, chl, a, N, load_sfx);
    kink_lst_fn = strcat(dn_lst(4), "/kink_lst.mat");
    if exist(kink_lst_fn, 'file') ~= 2
        fprintf("%s: '%s' does not exist.\n", datestr(datetime('now')), kink_lst_fn);
        return;
    end

    load(kink_lst_fn, "kink_lst");

    if(contains(load_str, 'u'))
        plot(x_lst, kink_lst(1, :), ...
            'LineStyle', ls_lst(1), ...
            'Color', color, ...
            'DisplayName', legstr_lst(1));
    end
    if(contains(load_str, 's'))
        plot(x_lst, kink_lst(2, :), ...
            'LineStyle', ls_lst(2), ...
            'Color', color, ...
            'DisplayName', legstr_lst(2));
    end
    if(contains(load_str, 'b'))
        plot_idx = 3 - (mode ~= 1);
        plot(x_lst, kink_lst(plot_idx, :), ...
            'LineStyle', ls_lst(plot_idx), ...
            'Color', color, ...
            'DisplayName', legstr_lst(plot_idx));
    end
    
    lg_obj = legend("show");
    set(lg_obj, 'Interpreter', 'latex');
end

function plot_kink_sc_rot(mode, load_str, color, x_lst)
%
% Plot kink of a single crack being rotated.
%

    chl = 1.0;
    
    ls_lst = [ "-" "-." "--" ];
    legstr_lst = ...
    [   ...
        "Uniaxial SC" ...
        "Shear SC" ...
        "Biaxial SC" ...
    ];
    if(mode ~= 1)
        chl = 2.0;
    end
    
    title("Single crack rotation");
    ylabel("Kinking angle $\theta_*$ [deg]");
    yticks(-180:20:180);
    xlabel("Rotation angle [deg]");
    
    if ~exist('x_lst','var')
        x_lst = 0:90;
    end        

    kk_sc = zeros(3, length(x_lst));
    
    for(i = 1:length(x_lst))
        if ~exist('x_lst','var')
            kink_mat = eval_kink_single_crack(chl, x_lst(i)*pi/180, 'ubs');
        else
            kink_mat = eval_kink_single_crack(chl, 0, 'ubs');
        end
        kk_sc(:, i) = kink_mat(1, :);
    end
    
    kk_sc(isnan(kk_sc)) = 0;
    
    if(contains(load_str, 'u'))
        plot(x_lst, kk_sc(1, :), ...
            'LineStyle', ls_lst(1), ...
            'Color', color, ...
            'DisplayName', legstr_lst(1));
    end
    if(contains(load_str, 's'))
        plot(x_lst, kk_sc(2, :), ...
            'LineStyle', ls_lst(2), ...
            'Color', color, ...
            'DisplayName', legstr_lst(2));
    end
    if(contains(load_str, 'b'))
        plot(x_lst, kk_sc(3, :), ...
            'LineStyle', ls_lst(3), ...
            'Color', color, ...
            'DisplayName', legstr_lst(3));
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
    if contains(load_str, 'b')
        [ k1, k2 ] = sif_tens_sc(1.0);
        
        sif_mat(1, 1:2, load_idx) = k1;
        sif_mat(1, 3:4, load_idx) = k2;
    end
    
    res = squeeze(eval_kink(sif_mat));    
end

function plot_stress(def_type, mode, a_coef, lstyle)

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
    clr_lst = ...
    [ ...
        "red" ...
        "green" ...
        "blue" ...
        "black" ...
    ];

    chl = 1;            % for input files
	a = a_coef*chl;
    N = 4;
    period_str = sprintf("$a = %.1fl$", a_coef);
    if(mode == 3)
        period_str = sprintf("$a = %.1fl_x$", a_coef);
    end
    case_lst = ...
    [   ...
        strcat("Uniaxial load. ", period_str) ...
        strcat("Shear load. ", period_str) ...
        strcat("Biaxial load. ", period_str) ...
    ];
    legstr_lst = [ "Crack 1"   "Crack 3" ...
                   "Crack 2" "Crack 4" ];

    if(mode ~= 1)
        case_lst(2) = case_lst(3);
    end    
    if(def_type == 1)
        x_lst = 0:90;
    else
        x_lst = 0 : 0.1*chl : (a/2);
        if mode == 3
            x_lst = 2*x_lst;
        end
    end
    
    if(mode == 1)
        title_pfx = "Square lattice. ";
        load_sfx = 'ubs';
    elseif(mode == 2)
        title_pfx = "Rectangular lattice. ";
        load_sfx = 'ub';
    else
        title_pfx = "Square lattice, $l_x = 2l_y$. ";
        load_sfx = 'ub';
    end
    
    if(mode ~= 1)
        chl_lst(1:2) = 2 * chl_lst(1:2);        % all horizontal are twice longer
    end
    
    dn_lst = kink_create_out_dirs(false, def_type, mode, chl, a, N, load_sfx);
    
    ss = zeros(size(x_lst));
    kk = zeros(size(x_lst));
    
    for(i_case = 1:strlength(load_sfx))
        figure('Name', strcat("Kinking angle at crack tips. ", title_pfx, case_lst(i_case)));        
        set(gcf,'color','w');
        box on;
        
        for(i_tip = 1:length(idx_lst))

            for(i = 1:length(x_lst))
                if(def_type == 1)
                    sif_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", x_lst(i)));
                    kink_fn = strcat(dn_lst(3), sprintf("/phi=%d.mat", x_lst(i)));
                else
                    sif_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", x_lst(i)));
                    kink_fn = strcat(dn_lst(3), sprintf("/yoff=%.3f.mat", x_lst(i)));
                end

                load(sif_fn, "sif_mat");
                load(kink_fn, "kink_mat");

                kk(i) = kink_mat(idx_lst(i_tip), tip_lst(i_tip), i_case);
                ss(i) = eval_stress(chl_lst(i_tip), sif_mat(idx_lst(i_tip), [tip_lst(i_tip), tip_lst(i_tip)+2], i_case), kk(i));
            end

            if(abs(kk(1)) > 170 && kk(1)*kk(2) < 0)
                kk(1) = (-1)*kk(1);
            end
            
            if(def_type == 2)               
                
                tmp(1, :) = x_lst;
                tmp(2, :) = kk;
                tmp(3, :) = ss;
                
                kk = spline(x_lst, kk, 0 : 0.01*chl : x_lst(end));
                ss = spline(x_lst, ss, 0 : 0.01*chl : x_lst(end));
                x_lst = 0 : 0.01*chl : x_lst(end);
            end
            
%             subplot(2, 1, 1);
%             hold on;
%             plot(x_lst, ss, ...
%                  'LineStyle', lstyle, ...
%                  'Color', clr_lst(i_tip), ...
%                  'DisplayName', legstr_lst(i_tip));
% 
%             subplot(2, 1, 2);
            hold on;
            plot(x_lst, kk, ...
                 'LineStyle', lstyle, ...
                 'Color', clr_lst(i_tip), ...
                 'DisplayName', legstr_lst(i_tip));
             
             if(def_type == 2)   
                x_lst = tmp(1, :);
                kk = tmp(2, :);
                ss = tmp(3, :);
            end
        end
        
%         subplot(2, 1, 1);
%         ylabel("Stress [Pa]");
%         yticks(0 : 0.1 : 1.6);
%         title(strcat("Relative stress at crack tips. ", title_pfx, case_lst(i_case)), 'Interpreter','latex');
%         lg_obj = legend("show");
%         set(lg_obj, 'Interpreter', 'latex');
%         if(def_type == 2)
%             xticks(x_lst);
%         end
%         
%         subplot(2, 1, 2);
        ylabel("Kinking angle $\theta_*$ [deg]");
        if(i_case < 3 - (strlength(load_sfx) == 2))
            yticks(-180 : 30 : 180);
        else
            yticks(-180 : 10 : 180);
        end
        if(def_type == 1)
            xlabel("Rotation angle $\alpha$ [deg]", 'Interpreter','latex');
            xticks(0:15:90);
        else
            xlabel("Relative uplift $\Delta / l$", 'Interpreter','latex');
            xticks(x_lst);
        end
        xlim([0 max(x_lst)]);
%         title(strcat("Kinking angle at crack tips. ", title_pfx, case_lst(i_case)), 'Interpreter','latex');
        lg_obj = legend("show");
        set(lg_obj, 'Interpreter', 'latex');
    end
end

function plot_err(form_type, def_type, mode, a_coef, lstyle)
%
% form_type = 1: plain ERR
% form_type = 2: ERR in kink direction

    function err = eval_ERR(sif)
    % sif: 
    %   2-vec with SIF I and II mode
    %%%
    
        k = sif(1);
        m = sif(2);
    
        if(form_type == 1)
            err = k^2 + m^2;
        else
            rt = sqrt(k^2 + 8*m^2);
            err = 128*m^4 * (k^2 + 6*m^2 + k*rt) / (16*m^2 + (k - rt)^2)^2;
        end
    end

    function s1 = eval_SIF1(sif, ka)
    % sif: 
    %   2-vec with SIF I and II mode
    % ka:
    %   kinking angle
    %%%
    
        k = sif(1);
        m = sif(2);
    
        if(form_type == 1)
            s1 = k;
        else
            c11 = (3*cos(ka/2) + cos(3*ka/2))/4;
            c12 =-3*(sin(ka/2) + sin(3*ka/2))/4;
            
            s1 = c11 * k + c12 * m;
        end
    end
    
    idx_lst = [ 4 2 9 10 ];                                                 % 4 2 9 10 | 2 9 10 1 13
    tip_lst = [ 1 2 2  1 ];                                                 % 1 2 2  1 | 1 1  2 1  2
    chl_lst = ones(length(idx_lst));
    clr_lst = ...
    [ ...
        "red" ...
        "green" ...
        "blue" ...
        "black" ...
    ];

    chl = 1;            % for input files
	a = a_coef*chl;
    N = 4;
    period_str = sprintf("$a = %.1fl$", a_coef);
    if(mode == 3)
        period_str = sprintf("$a = %.1fl_x$", a_coef);
    end
    case_lst = ...
    [   ...
        strcat("Uniaxial load. ", period_str) ...
        strcat("Shear load. ", period_str) ...
        strcat("Biaxial load. ", period_str) ...
    ];
    legstr_lst = [ "Crack 1"   "Crack 3" ...
                   "Crack 2" "Crack 4" ];

    if(mode ~= 1)
        case_lst(2) = case_lst(3);
    end    
    if(def_type == 1)
        x_lst = 0:90;
    else
        x_lst = 0 : 0.1*chl : (a/2);
        if mode == 3
            x_lst = 2*x_lst;
        end
    end
    
    if(form_type == 1)
        title_form_str = "";
    else
        title_form_str = " (at kink tip)";
    end
    
    if(mode == 1)
        title_pfx = "Square lattice. ";
        load_sfx = 'ubs';
    elseif(mode == 2)
        title_pfx = "Rectangular lattice. ";
        load_sfx = 'ub';
    else
        title_pfx = "Square lattice, $l_x = 2l_y$. ";
        load_sfx = 'ub';
    end
    
    dn_lst = kink_create_out_dirs(false, def_type, mode, chl, a, N, load_sfx);
    
    err = zeros(size(x_lst));
    k1 = zeros(size(x_lst));
    
    for(i_case = 1:strlength(load_sfx))
        fig_err = figure('Name', sprintf("Rel. ERR%s. %s%s", title_form_str, title_pfx, case_lst(i_case)));        
        set(gcf,'color','w');
        box on;
        
        fig_sif = figure('Name', sprintf("Rel. SIF mode I%s. %s%s", title_form_str, title_pfx, case_lst(i_case)));      
        set(gcf,'color','w');
        box on;
        
        for(i_tip = 1:length(idx_lst))

            for(i = 1:length(x_lst))
                if(def_type == 1)
                    sif_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", x_lst(i)));
                    kink_fn = strcat(dn_lst(3), sprintf("/phi=%d.mat", x_lst(i)));
                else
                    sif_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", x_lst(i)));
                    kink_fn = strcat(dn_lst(3), sprintf("/yoff=%.3f.mat", x_lst(i)));
                end
                
                load(sif_fn, "sif_mat");
                load(kink_fn, "kink_mat");
                
                err(i) = eval_ERR(sif_mat(idx_lst(i_tip), [tip_lst(i_tip), tip_lst(i_tip)+2], i_case)) / (pi*chl*(1 + (mode>1)));
                k1(i) = eval_SIF1(sif_mat(idx_lst(i_tip), [tip_lst(i_tip), tip_lst(i_tip)+2], i_case), ...
                                  kink_mat(idx_lst(i_tip), tip_lst(i_tip), i_case) * pi/180) / sqrt(pi*chl*(1 + (mode>1)));
            end
            
            if(def_type == 2)               
                
                tmp(1, :) = x_lst;
                tmp(2, :) = err;
                tmp(3, :) = k1;
                
                err = spline(x_lst, err, 0 : 0.01*chl : x_lst(end));
                k1 = spline(x_lst, k1, 0 : 0.01*chl : x_lst(end));
                x_lst = 0 : 0.01*chl : x_lst(end);
            end
            
            [maxerr, maxerr_idx] = max(err);
            [maxk1, maxk1_idx] = max(k1);
            
            figure(fig_err);
            hold on;
            plot(x_lst, err, ...
                 'LineStyle', lstyle, ...
                 'Color', clr_lst(i_tip), ...
                 'DisplayName', strcat(legstr_lst(i_tip))); %, sprintf(". Max at %.3f", x_lst(maxerr_idx))));

            figure(fig_sif);
            hold on;
            plot(x_lst, k1, ...
                 'LineStyle', lstyle, ...
                 'Color', clr_lst(i_tip), ...
                 'DisplayName', strcat(legstr_lst(i_tip))); %, sprintf(". Max at %.3f", x_lst(maxk1_idx))));
             
            if(def_type == 2)   
                x_lst = tmp(1, :);
                err = tmp(2, :);
                k1 = tmp(3, :);
            end
        end
        
        figure(fig_err);
        ylabel("$G(\theta_*) / G^0$", 'Interpreter', 'latex');
%         if(i_case < 3 - (strlength(load_sfx) == 2))
%             yticks(-180 : 30 : 180);
%         else
%             yticks(-180 : 10 : 180);
%         end
        if(def_type == 1)
            xlabel("Rotation angle $\alpha$ [deg]", 'Interpreter','latex');
            xticks(0:15:90);
        else
            xlabel("Relative uplift $\Delta / l$", 'Interpreter','latex');
            xticks(x_lst);
        end
        xlim([0 max(x_lst)]);
        
        lg_obj = legend("show");
        set(lg_obj, 'Interpreter', 'latex');
        
        figure(fig_sif);
        ylabel("$K_I / K_I^0$", 'Interpreter', 'latex');
        if(def_type == 1)
            xlabel("Rotation angle $\alpha$ [deg]", 'Interpreter','latex');
            xticks(0:15:90);
        else
            xlabel("Relative uplift $\Delta / l$", 'Interpreter','latex');
            xticks(x_lst);
        end
        xlim([0 max(x_lst)]);
        
        lg_obj = legend("show");
        set(lg_obj, 'Interpreter', 'latex');
    end
end

function plot_tstress(mode, a_coef, color)

    chl = 1;
    a = a_coef*chl;
    N = 4;
    load_sfx = "ub";
    
    ls_lst = [ "-" "--" ];
    period_str = sprintf("$a = %.1fl$", a_coef);
    if(mode == 3)
        period_str = sprintf("$a = %.1fl_x$", a_coef);
    end

    ylabel("Kinking angle $\theta_*$ [deg]");
%     yticks(-180:20:180);
%     if(mode ~= 1)
%         yticks(-180:10:180);
%     end
    
    xlabel("$\lambda = \sigma^\infty_{xx} / \sigma^\infty_{yy}$", 'Interpreter', 'latex');
    x_lst = -1 : 0.1 : 1;
    
    if mode == 1
        leg_pfx = "Square";
    elseif mode == 2
        leg_pfx = "Rect";
    else
        leg_pfx = "Square, $l_x = 2l_y$";
    end
 
    dn_lst = kink_create_out_dirs(false, 3, mode, chl, a, N, load_sfx);
    kink_lst_fn = strcat(dn_lst(4), "/kink_lst.mat");
    kink_new_lst_fn = strcat(dn_lst(4), "/kink_new_lst.mat");
    if exist(kink_lst_fn, 'file') ~= 2
        fprintf("%s: '%s' does not exist.\n", datestr(datetime('now')), kink_lst_fn);
        return;
    end

    load(kink_lst_fn, "kink_lst");
    load(kink_new_lst_fn, "kink_new_lst");
    kink_new_lst = kink_new_lst * 180/pi;
    
    leg1 = sprintf("Old. %s. %s", leg_pfx, period_str);
    leg2 = sprintf("New. %s. %s", leg_pfx, period_str);
    
    plot(x_lst, kink_lst, 'LineStyle', ls_lst(1), 'Color', color, 'DisplayName', leg1);
    plot(x_lst, kink_new_lst, 'LineStyle', ls_lst(2), 'Color', color, 'DisplayName', leg2);
    
    lg_obj = legend("show");
    set(lg_obj, 'Interpreter', 'latex');
end

function plot_sif3(def_type, mode, a_coef, lstyle)
    
    idx_lst = [ 4 2 9 10 ];                                                 % 4 2 9 10 | 2 9 10 1 13
    tip_lst = [ 1 2 2  1 ];                                                 % 1 2 2  1 | 1 1  2 1  2
    chl_lst = ones(length(idx_lst));
    clr_lst = ...
    [ ...
        "red" ...
        "green" ...
        "blue" ...
        "black" ...
    ];

    chl = 1;            % for input files
	a = a_coef*chl;
    N = 4;
    load_sfx = 'y';
    period_str = sprintf("$a = %.1fl$", a_coef);
    if(mode == 3)
        period_str = sprintf("$a = %.1fl_x$", a_coef);
    end
    
    legstr_lst = [ "Crack 1"   "Crack 3" ...
                   "Crack 2" "Crack 4" ];

    if(def_type == 1)
        x_lst = 0:90;
    else
        x_lst = 0 : 0.1*chl : (a/2);
        if mode == 3
            x_lst = 2*x_lst;
        end
    end
    
    if(mode == 1)
        title_pfx = "Square lattice. ";
    elseif(mode == 2)
        title_pfx = "Rectangular lattice. ";
    else
        title_pfx = "Square lattice, $l_x = 2l_y$. ";
    end
    
    dn_lst = sif3_create_out_dirs(false, def_type, mode, chl, a, N, load_sfx);
    
    k3 = zeros(size(x_lst));
    
    for(i_tip = 1:length(idx_lst))

        for(i = 1:length(x_lst))
            if(def_type == 1)
                sif_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", x_lst(i)));
            else
                sif_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", x_lst(i)));
            end

            load(sif_fn, "sif_mat");

            k3(i) = sif_mat(idx_lst(i_tip), tip_lst(i_tip)) / sqrt(pi*chl*(1 + (mode==3)));
        end

        hold on;
        [maxk3, max_idx] = max(k3);
        plot(x_lst, k3, ...
             'LineStyle', lstyle, ...
             'Color', clr_lst(i_tip), ...
             'DisplayName', strcat(legstr_lst(i_tip))); %, sprintf(". Max, Argmax = %.3f, %.1f", maxk3, x_lst(max_idx))));
    end

%     if(i_case < 3 - (strlength(load_sfx) == 2))
%         yticks(-180 : 30 : 180);
%     else
%         yticks(-180 : 10 : 180);
%     end
    ylabel("$K_{III} / K_{III}^0$", 'Interpreter', 'latex');
    if(def_type == 1)
        xlabel("Rotation angle $\alpha$ [deg]", 'Interpreter','latex');
    else
        xlabel("Relative uplift $\Delta / l$", 'Interpreter','latex');
        xticks(x_lst);
    end
    xlim([0 max(x_lst)]);
%     title(strcat("SIF mode III. Antiplane ZY. ", title_pfx, period_str), 'Interpreter','latex');
    lg_obj = legend("show");
    set(lg_obj, 'Interpreter', 'latex');
    set(gcf,'color','w');
    box on;
end







































