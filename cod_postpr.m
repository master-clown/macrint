function cod_postpr()
    
%     plot_cod(1, 1, 2.5);
%     plot_cod(1, 1, 3);
%     plot_cod(1, 1, 4);
%     plot_cod(1, 2, 3);
%     plot_cod(1, 3, 3);
% 
%     plot_cod(2, 1, 2.5);
%     plot_cod(2, 1, 3);
%     plot_cod(2, 1, 4);
%     plot_cod(2, 2, 3);
%     plot_cod(2, 3, 3);    

    N_cell_h_cr = 4;
    chly = 1;
    chlx = chly / sqrt(N_cell_h_cr);
    ay = 2.4 * chly;
   %ay = 4.8 * chly;
    ax = 2 * ay;
   %ax = 4 * ay;
   %ax = 1 * ay;
    
    eval_edge_effect(ax, ay, chlx, chly, N_cell_h_cr);
    
end

function plot_cod(def_type, mode, a_coef)
    
    idx_lst = [ 4 9 2 10 5 3 ];                                            
    chl_lst = ones(length(idx_lst));
    sty_lst = [ "-", "--" ];                % components
    clr_lst = ...                           % load types
    [ ...
        "red" ...
        "green" ...
        "blue" ...
    ];

    chl = 1;            % for input files
	a = a_coef*chl;
    N = 4;
    period_str = sprintf("$a = %.1fl$", a_coef);
    if(mode == 3)
        period_str = sprintf("$a = %.1fl_x$", a_coef);
    end
    
    comp_str_lst = [ "Normal comp", "Tangent comp" ];
    loty_str_lst = ...
    [   ...
        "Uniax" ...
        "Shear" ...
        "Biax" ...        
    ];    
    if(mode ~= 1)
        loty_str_lst(2) = loty_str_lst(3);
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
        fig_title = strcat("COD. Square lattice. ", period_str);
        load_sfx = 'ubs';
    elseif(mode == 2)
        fig_title = strcat("COD. Rectangular lattice. ", period_str);
        load_sfx = 'ub';
    else
        fig_title = strcat("COD. Square lattice, $l_x = 2l_y$. ", period_str);
        load_sfx = 'ub';
    end
    
    if(mode ~= 1)
        chl_lst([1 3 5 6]) = 2 * chl_lst([1 3 5 6]);        % all horizontal are twice longer
    end
    
    dn_lst = kink_create_out_dirs(false, def_type, mode, chl, a, N, load_sfx);
    
    figure('Name', fig_title);
    set(gcf,'color','w');
    box on;
    
    cod_m = zeros(length(idx_lst), 2, strlength(load_sfx), length(x_lst));     % [...] X [normal and tang comps] X [ load types] X [...]
    
    for(i  = 1:length(x_lst))
        if(def_type == 1)
            avg_tr_fn = strcat(dn_lst(2), sprintf("/phi=%d.mat", x_lst(i)));
        else
            avg_tr_fn = strcat(dn_lst(2), sprintf("/yoff=%.3f.mat", x_lst(i)));
        end

        load(avg_tr_fn, "avg_tr");
        
        for(i_cr = 1:length(idx_lst))
            cod_m(i_cr, :, :, i) = pi*chl_lst(i_cr) * avg_tr(idx_lst(i_cr), :, :);
        end
    end

    for(i_cr = 1:length(idx_lst))

        subplot(2, 3, i_cr);
        hold on;
        title(sprintf("Crack %d", i_cr));
        
        for(i_comp = 1:2)
        for(i_loty = 1:strlength(load_sfx))
            
            legstr = sprintf("%s. %s", comp_str_lst(i_comp), loty_str_lst(i_loty));
            
            plot(x_lst, squeeze(cod_m(i_cr, i_comp, i_loty, :)), ...
                 'LineStyle', sty_lst(i_comp), ...
                 'Color', clr_lst(i_loty), ...
                 'DisplayName', legstr);
        end
        end
        
        if(def_type == 1)
            %xlabel("Rotation angle $\alpha$ [deg]", 'Interpreter','latex');
            xticks(0:15:90);
        else
            %xlabel("Relative uplift $\Delta / l$", 'Interpreter','latex');
            xticks(x_lst);
        end
        xlim([0 max(x_lst)]);
        
        lg_obj = legend("show");
        set(lg_obj, 'Interpreter', 'latex');
    end
end

function eval_edge_effect(ax, ay, chlx, chly, N_cell_h_cr)

    load_sfx = 'ubs';
    
    err = zeros(2, 3);
    
    cr_idx = 30;
    for N = 3:4

        %измени имя папки откуда читаются результаты
        dn_lst = cod_create_out_dirs(0, ax, ay, chlx, chly, N_cell_h_cr, N, load_sfx);
        sif_fn = strcat(dn_lst(1), "/", "sif_mat.mat");
        avg_fn = strcat(dn_lst(1), "/", "avg_tr.mat");
        
        load(sif_fn, "sif_mat");
        load(avg_fn, "avg_tr");
        
        err(N-2, :) = sif_mat(cr_idx, 1, :).^2 + sif_mat(cr_idx, 3, :).^2;
    end
    
    ef = (err(2,:) - err(1,:)) ./ err(1,:);
    
    %тут перепутаны биакс и шир
    fprintf("EDGE EFFECT 16->20:\nUniax: %f\nBiax: %f\nShear: %f", ef(1), ef(3), ef(2));
end


























