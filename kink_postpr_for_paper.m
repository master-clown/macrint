function kink_postpr_for_paper()
    
    %% Common variables
    idx_lst = [ 4 2 9 10 ];                                                 % 4 2 9 10 | 2 9 10 1 13 -- 1 3 2 4
    tip_lst = [ 1 2 2  1 ];                                                 % 1 2 2  1 | 1 1  2 1  2
    cr_ID_lst = [ 1 3 2 4 ];
    legstr_lst = [ "Crack 1" "Crack 3" "Crack 2" "Crack 4" ];
    clr_lst = [ "green" "red" "blue", "black" ];
    ls_p = "-";
    ls_r = "--";
    ls_l = "none";    
    
    chl = 1;
    a = 2.5*chl;
    N = 4;
    load_sfx = 'ubs';
    loadstr_lst = [ "Uniaxial vertical load" "Shear load" "Hydrostatic load" ];
    ang_v = 0:90;
    off_v = 0 : 0.1*chl : (a/2);
    
    defl_xlab_lst = [ "Rotation angle $\alpha$ [deg]" "Relative uplift $\Delta / l$" ];
    defl_xtick_lst = { 0:15:90, off_v };
    
    set(0,'defaultAxesFontSize',12)
    
    %% Plots for collinear mode 4
    kink_postpr_for_paper_mode4();
    return;
    
%     %% Figure 4 and 6: SIF I, II, ERR and kinking angle, a = 2.5L
%     load_idx = 1;
%     x_v = ang_v;
%     for(i_defl_type = 1:2)
%         
%         % 'p' - by C.-Rice and Panasuk,'r' - by C.-Rice and maxERR, 'l' - by Leblond
%         sp_vals = zeros(4, length(x_v), 2);            % SIF I and II
%         ep_vals = zeros(4, length(x_v));               % ERR        
%         kp_vals = zeros(4, length(x_v));               % Kinking
%         sr_vals = zeros(4, length(x_v), 2);
%         er_vals = zeros(4, length(x_v));
%         kr_vals = zeros(4, length(x_v));                        
%         sl_vals = zeros(4, length(x_v), 2);
%         el_vals = zeros(4, length(x_v));
%         kl_vals = zeros(4, length(x_v));                          
%         for(i = 1:length(x_v))
% 
%             dn_lst = kink_create_out_dirs(false, i_defl_type, 1, chl, a, N, load_sfx);
%             if(i_defl_type == 1)
%                 sif_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", x_v(i)));
%             else
%                 sif_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", x_v(i)));
%             end
%             load(sif_fn, "sif_mat");
%             
%             [sp_m,kp_m] = eval_kink_panasuk(sif_mat);                           % Panasuk
%             [sr_m,kr_m] = eval_kink_maxerr(sif_mat);                            % Rice
%             [sl_m,kl_m] = eval_kink_lebl(sif_mat);
% 
%             for(i_cr = 1:4)
%                 kp_vals(i_cr, i) = kp_m(idx_lst(i_cr), tip_lst(i_cr), load_idx);
%                 kr_vals(i_cr, i) = kr_m(idx_lst(i_cr), tip_lst(i_cr), load_idx);
%                 kl_vals(i_cr, i) = kl_m(idx_lst(i_cr), tip_lst(i_cr), load_idx);
%                 
%                 sp_vals(i_cr, i, 1) = sp_m(idx_lst(i_cr), tip_lst(i_cr), load_idx) / sqrt(pi*chl);
%                 sr_vals(i_cr, i, 1) = sr_m(idx_lst(i_cr), tip_lst(i_cr), load_idx) / sqrt(pi*chl);
%                 sl_vals(i_cr, i, 1) = sl_m(idx_lst(i_cr), tip_lst(i_cr), load_idx) / sqrt(pi*chl);
%                 sp_vals(i_cr, i, 2) = sp_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx) / sqrt(pi*chl);
%                 sr_vals(i_cr, i, 2) = sr_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx) / sqrt(pi*chl);
%                 sl_vals(i_cr, i, 2) = sl_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx) / sqrt(pi*chl);
%                 
%                 ep_vals(i_cr, i) = (sp_m(idx_lst(i_cr), tip_lst(i_cr), load_idx)^2 + ...
%                                     sp_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx)^2) / (pi*chl);   
%                 er_vals(i_cr, i) = (sr_m(idx_lst(i_cr), tip_lst(i_cr), load_idx)^2 + ...
%                                     sr_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx)^2) / (pi*chl);   
%                 el_vals(i_cr, i) = (sl_m(idx_lst(i_cr), tip_lst(i_cr), load_idx)^2 + ...
%                                     sl_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx)^2) / (pi*chl);                
%             end
%             
%         end
%         
%         %% %% Plot ERR
%         title_lst = [ sprintf("Fig %da. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ...
%                       sprintf("Fig %db. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ];
%         fn_lst = [ sprintf("plots/Fig_%da", 7 - load_idx), sprintf("plots/Fig_%db", 7 - load_idx) ];
%         
%         figure('Name', strcat(title_lst(i_defl_type), loadstr_lst(load_idx)), 'Color', 'w');
%         
%         hold on;
%         box on;
%         
%         cr_to_plot = 4 - 2*(load_idx == 1);
%         for(i_cr = 1:cr_to_plot)
%             vals_p = ep_vals(i_cr, :);
%             vals_r = er_vals(i_cr, :);
%             vals_l = el_vals(i_cr, :);
%             if(i_defl_type == 2)
%                 x_v = 0 : 0.01*chl : x_v(end);
%                 vals_p = spline(off_v, ep_vals(i_cr, :), x_v);
%                 vals_r = spline(off_v, er_vals(i_cr, :), x_v);
%                 vals_l = spline(off_v, el_vals(i_cr, :), x_v);
%             end
%             
%             plot(x_v, vals_p, ...
%                  'LineStyle', ls_p, ...
%                  'Color', clr_lst(i_cr), ...
%                  'LineWidth', 1.5);
%             plot(x_v, vals_r, ...
%                  'LineStyle', ls_r, ...
%                  'Color', clr_lst(i_cr));
%             
%             if(i_defl_type == 1)
%                 idx_lebl = 1:4:91;            % every 5 degrees
%             else
%                 idx_lebl = 1:10:121;     
%             end
%             
%             plot(x_v(idx_lebl), vals_l(idx_lebl), ...
%                  'LineStyle', ls_l, ...
%                  'Color', clr_lst(i_cr), ...
%                  'Marker', 'o');
%              
%             Xl = xlim();
%             Yl = ylim();
%             xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
%             ybeg = (vals_p(50) - Yl(1)) / (Yl(2) - Yl(1));
%             
%             annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)));
%             
%         end
%         
%         xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
%         xticks(defl_xtick_lst{i_defl_type});
%         xlim([0 max(x_v)]);
%         ylabel("$G(\theta_*) / G^0$", 'Interpreter', 'latex');
%         %yticks(-90:20:90);
%         
% %         lg_obj = legend("show");
% %         set(lg_obj, 'Interpreter', 'latex');
% 
%         savefig(fn_lst(i_defl_type));
%         
%         %% %% Plot SIF I
%         title_lst = [ sprintf("Fig %dc. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ...
%                       sprintf("Fig %dd. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ];
%         fn_lst = [ sprintf("plots/Fig_%dc", 7 - load_idx), sprintf("plots/Fig_%dd", 7 - load_idx) ];
%         
%         figure('Name', strcat(title_lst(i_defl_type), loadstr_lst(load_idx)), 'Color', 'w');
%         hold on;
%         box on;
%         
%         cr_to_plot = 4 - 2*(load_idx == 1);
%         for(i_cr = 1:cr_to_plot)
%             vals_p = sp_vals(i_cr, :, 1);
%             vals_r = sr_vals(i_cr, :, 1);
%             vals_l = sl_vals(i_cr, :, 1);
%             if(i_defl_type == 2)
%                 x_v = 0 : 0.01*chl : x_v(end);
%                 vals_p = spline(off_v, sp_vals(i_cr, :, 1), x_v);
%                 vals_r = spline(off_v, sr_vals(i_cr, :, 1), x_v);
%                 vals_l = spline(off_v, sl_vals(i_cr, :, 1), x_v);
%             end
%             
%             plot(x_v, vals_p, ...
%                  'LineStyle', ls_p, ...
%                  'Color', clr_lst(i_cr), ...
%                  'LineWidth', 1.5);
%             plot(x_v, vals_r, ...
%                  'LineStyle', ls_r, ...
%                  'Color', clr_lst(i_cr));
%             
%             if(i_defl_type == 1)
%                 idx_lebl = 1:4:91;            % every 5 degrees
%             else
%                 idx_lebl = 1:10:121;     
%             end
%             
%             plot(x_v(idx_lebl), vals_l(idx_lebl), ...
%                  'LineStyle', ls_l, ...
%                  'Color', clr_lst(i_cr), ...
%                  'Marker', 'o');
%              
%             Xl = xlim();
%             Yl = ylim();
%             xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
%             ybeg = (vals_p(50) - Yl(1)) / (Yl(2) - Yl(1));
%             
%             annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         end
%         
%         xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
%         xticks(defl_xtick_lst{i_defl_type});
%         xlim([0 max(x_v)]);
%         ylabel("$\hat{K}_I / K_I^0$", 'Interpreter', 'latex');
%         %yticks(-90:20:90);
%         
%         savefig(fn_lst(i_defl_type));
%         
% %         lg_obj = legend("show");
% %         set(lg_obj, 'Interpreter', 'latex');
%         
% %         %% Plot SIF II
% %         title_lst = [ "Fig 6g. Square, a = 2.5. " ...
% %                       "Fig 6h. Square, a = 2.5. " ];
% %         figure('Name', strcat(title_lst(i_defl_type), loadstr_lst(load_idx)), 'Color', 'w');
% %         hold on;
% %         
% %         for(i_cr = 1:4)
% %             vals1 = sr_vals(i_cr, :, 2);
% %             vals2 = sl_vals(i_cr, :, 2);
% %             if(i_defl_type == 2)
% %                 x_v = 0 : 0.01*chl : x_v(end);
% %                 vals1 = spline(off_v, sr_vals(i_cr, :, 2), x_v);
% %                 vals2 = spline(off_v, sl_vals(i_cr, :, 2), x_v);
% %             end
% %             
% %             plot(x_v, vals1, ...
% %                  'LineStyle', ls_lst(i_cr), ...
% %                  'Color', clr_lst(i_cr), ...
% %                  'DisplayName', strcat(legstr_lst(i_cr), ". Rice"));
% %             plot(x_v, vals2, ...
% %                  'LineStyle', ls_lst(i_cr), ...
% %                  'Color', clr_lst(i_cr), ...
% %                  'LineWidth', 2, ...
% %                  'DisplayName', strcat(legstr_lst(i_cr), ". Leblond"));
% %         end
% %         
% %         xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
% %         xticks(defl_xtick_lst{i_defl_type});
% %         xlim([0 max(x_v)]);
% %         ylabel("$\hat{K}_{II} / K_{I}^0$", 'Interpreter', 'latex');
% %         %yticks(-90:20:90);
% %         
% %         lg_obj = legend("show");
% %         set(lg_obj, 'Interpreter', 'latex');
%         
%         %% Plot kinking angle
%         title_lst = [ sprintf("Fig %de. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ...
%                       sprintf("Fig %df. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ];
%         fn_lst = [ sprintf("plots/Fig_%de", 7 - load_idx), sprintf("plots/Fig_%df", 7 - load_idx) ];
%         
%         figure('Name', strcat(title_lst(i_defl_type), loadstr_lst(load_idx)), 'Color', 'w');
%         hold on;
%         box on;
%         
%         cr_to_plot = 4 - 2*(load_idx == 1);
%         for(i_cr = 1:cr_to_plot)
%             vals_p = kp_vals(i_cr, :);
%             vals_r = kr_vals(i_cr, :);
%             vals_l = kl_vals(i_cr, :);
%             if(i_defl_type == 2)
%                 x_v = 0 : 0.01*chl : x_v(end);
%                 vals_p = spline(off_v, kp_vals(i_cr, :), x_v);
%                 vals_r = spline(off_v, kr_vals(i_cr, :), x_v);
%                 vals_l = spline(off_v, kl_vals(i_cr, :), x_v);
%             end
%             
%             plot(x_v, vals_p, ...
%                  'LineStyle', ls_p, ...
%                  'Color', clr_lst(i_cr), ...
%                  'LineWidth', 1.5);
%             plot(x_v, vals_r, ...
%                  'LineStyle', ls_r, ...
%                  'Color', clr_lst(i_cr));
%             
%             if(i_defl_type == 1)
%                 idx_lebl = 1:4:91;            % every 5 degrees
%             else
%                 idx_lebl = 1:10:121;     
%             end
%             
%             plot(x_v(idx_lebl), vals_l(idx_lebl), ...
%                  'LineStyle', ls_l, ...
%                  'Color', clr_lst(i_cr), ...
%                  'Marker', 'o');
%              
%             Xl = xlim();
%             Yl = ylim();
%             xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
%             ybeg = (vals_p(50) - Yl(1)) / (Yl(2) - Yl(1));
%             
%             annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         end
%         
%         xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
%         xticks(defl_xtick_lst{i_defl_type});
%         xlim([0 max(x_v)]);
%         ylabel("Kinking angle $\theta_*$ [deg]", 'Interpreter', 'latex');
%         yticks(-90:5:90);
%         
% %         Leg1 = plot(nan, nan, 'k-');
% %         Leg2 = plot(nan, nan, 'k--');
% %         Leg3 = plot(nan, nan, 'LineStyle', ls_l, 'Marker', 'o', 'Color', 'black');
% %         legend([Leg1 Leg2 Leg3], 'One', 'Two', 'Three' );
%         
%         savefig(fn_lst(i_defl_type));
% 
%         x_v = off_v;
%     end
%     
%     
%     %% Figure 3 or 5: ERR, SIF I and II for parent cracks
% %     load_idx = 1;
%     x_v = ang_v;
%     for(i_defl_type = 1:2)
%         
%         s1_vals = zeros(4, length(x_v));            % SIF I                     
%         s2_vals = zeros(4, length(x_v));            % SIF II
%         err_vals= zeros(4, length(x_v));
%         for(i = 1:length(x_v))
% 
%             dn_lst = kink_create_out_dirs(false, i_defl_type, 1, chl, a, N, load_sfx);
%             if(i_defl_type == 1)
%                 sif_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", x_v(i)));
%             else
%                 sif_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", x_v(i)));
%             end
%             load(sif_fn, "sif_mat");
% 
%             for(i_cr = 1:4)
%                 s1_vals(i_cr, i) = sif_mat(idx_lst(i_cr), tip_lst(i_cr), load_idx) / sqrt(pi*chl);
%                 s2_vals(i_cr, i) = sif_mat(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx) / sqrt(pi*chl);
%                 err_vals(i_cr, i) = (s1_vals(i_cr, i)^2 + s2_vals(i_cr, i)^2);
%             end
%             
%         end
%         
%         %% %% Plot SIF I
%         title_lst = [ sprintf("Fig %da. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ...
%                       sprintf("Fig %db. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ];
%         fn_lst = [ sprintf("plots/Fig_%da", 6 - load_idx), sprintf("plots/Fig_%db", 6 - load_idx) ];
%         
%         figure('Name', title_lst(i_defl_type), 'Color', 'w');
%         hold on;
%         box on;
%         
%         for(i_cr = 1:4)
%             vals1 = s1_vals(i_cr, :);
%             if(i_defl_type == 2)
%                 x_v = 0 : 0.01*chl : x_v(end);
%                 vals1 = spline(off_v, s1_vals(i_cr, :), x_v);
%             end
%             
%             plot(x_v, vals1, ...
%                  'LineStyle', '-', ...
%                  'Color', clr_lst(i_cr));
%             
%             Xl = xlim();
%             Yl = ylim();
%             xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
%             ybeg = (vals1(50) - Yl(1)) / (Yl(2) - Yl(1));
%             
%             annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         end
%         
%         xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
%         xticks(defl_xtick_lst{i_defl_type});
%         xlim([0 max(x_v)]);
%         ylabel("$K_I / K_I^0$", 'Interpreter', 'latex');
%         %yticks(-90:20:90);
%         
% %         lg_obj = legend("show");
% %         set(lg_obj, 'Interpreter', 'latex');
% 
%         savefig(fn_lst(i_defl_type));
%         
%         %% Plot SIF II
%         title_lst = [ sprintf("Fig %dc. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ...
%                       sprintf("Fig %dd. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ];
%         fn_lst = [ sprintf("plots/Fig_%dc", 6 - load_idx), sprintf("plots/Fig_%dd", 6 - load_idx) ];
%         
%         figure('Name', title_lst(i_defl_type), 'Color', 'w');
%         hold on;
%         box on;
%         
%         for(i_cr = 1:4)
%             vals1 = s2_vals(i_cr, :);
%             if(i_defl_type == 2)
%                 x_v = 0 : 0.01*chl : x_v(end);
%                 vals1 = spline(off_v, s2_vals(i_cr, :), x_v);
%             end
%             
%             plot(x_v, vals1, ...
%                  'LineStyle', '-', ...
%                  'Color', clr_lst(i_cr));
%              
%             Xl = xlim();
%             Yl = ylim();
%             xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
%             ybeg = (vals1(50) - Yl(1)) / (Yl(2) - Yl(1));
%             
%             annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         end
%         
%         xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
%         xticks(defl_xtick_lst{i_defl_type});
%         xlim([0 max(x_v)]);
%         ylabel("$K_{II} / K_{I}^0$", 'Interpreter', 'latex');
%         %yticks(-90:20:90);
%         
% %         lg_obj = legend("show");
% %         set(lg_obj, 'Interpreter', 'latex');
% 
%         savefig(fn_lst(i_defl_type));
%         
%         %% Plot ERR
%         title_lst = [ sprintf("Fig %de. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ...
%                       sprintf("Fig %df. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ];
%         fn_lst = [ sprintf("plots/Fig_%de", 6 - load_idx), sprintf("plots/Fig_%df", 6 - load_idx) ];
%         
%         figure('Name', title_lst(i_defl_type), 'Color', 'w');
%         hold on;
%         box on;
%         
%         for(i_cr = 1:4)
%             vals1 = err_vals(i_cr, :);
%             if(i_defl_type == 2)
%                 x_v = 0 : 0.01*chl : x_v(end);
%                 vals1 = spline(off_v, err_vals(i_cr, :), x_v);
%             end
%             
%             plot(x_v, vals1, ...
%                  'LineStyle', '-', ...
%                  'Color', clr_lst(i_cr));
%             
%             Xl = xlim();
%             Yl = ylim();
%             xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
%             ybeg = (vals1(50) - Yl(1)) / (Yl(2) - Yl(1));
%             
%             annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         end
%         
%         xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
%         xticks(defl_xtick_lst{i_defl_type});
%         xlim([0 max(x_v)]);
%         ylabel("$G / G^0$", 'Interpreter', 'latex');
%         %yticks(-90:20:90);
%         
% %         lg_obj = legend("show");
% %         set(lg_obj, 'Interpreter', 'latex');
%         
%         savefig(fn_lst(i_defl_type));
% 
%         x_v = off_v;
%     end
    
    %% Period influence analysis
    load_idx = 1;
    x_v = ang_v;
    
    period_lst = [ 2.5, 3, 4, 3, 3 ]; 
    geom_type_lst = [ 1, 1, 1, 2, 3 ];
    load_sfx_lst = [ "ubs" "ubs" "ubs" "ub" "ub" ];
    ls_lst = [ "--" "-" ":" "--" ":"];
    thc_lst = [ 1.5 1.5 1.5 3 ];
    loadstr_lst = [ "Uniaxial vertical load" "Hydrostatic load" ];
    x_off_v_lst = { 0 : 0.1*chl :  period_lst(1)/2, ...
                    0 : 0.1*chl :  period_lst(2)/2, ...
                    0 : 0.1*chl :  period_lst(3)/2, ...
                    0 : 0.1*chl :  period_lst(4)/2, ...
                    0 : 0.2*chl :  period_lst(5) };
    load_idx_lst = [ 1 1 1 1 1;
                     3 3 3 2 2 ];
       
    
    for(i_defl_type = 1:2)
        
        if(i_defl_type == 1)
            x_siz = length(x_v);
        else
            x_siz = 21;
        end
        
        % 'p' - by C.-Rice and Panasuk,'r' - by C.-Rice and maxERR, 'l' - by Leblond
        % there are 5 periods
        ep_vals = zeros(4, x_siz, 5);               % ERR        
        kp_vals = zeros(4, x_siz, 5);               % Kinking
        er_vals = zeros(4, x_siz, 5);
        kr_vals = zeros(4, x_siz, 5);         
        el_vals = zeros(4, x_siz, 5);
        kl_vals = zeros(4, x_siz, 5);                          
        for(i_period = 1:5)

            load_ind = load_idx_lst(load_idx, i_period);
            x_off_v = x_off_v_lst{i_period};
            x_off_len = length(x_off_v);
            dn_lst = kink_create_out_dirs(false, i_defl_type, geom_type_lst(i_period), chl, period_lst(i_period), N, load_sfx_lst(i_period));
            
            for(i = 1:x_siz)                
                if(i_defl_type == 1)
                    sif_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", x_v(i)));
                else
                    if i > x_off_len
                        break;
                    end
                    sif_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", x_off_v(i)));
                end
                load(sif_fn, "sif_mat");

                [sp_m,kp_m] = eval_kink_panasuk(sif_mat);                           % Panasuk
                [sr_m,kr_m] = eval_kink_maxerr(sif_mat);                            % Rice
                [sl_m,kl_m] = eval_kink_lebl(sif_mat);

                coef = [ 1 1 1 1 ];
                if i_period == 4 || i_period == 5
                    coef = [ 2 2 1 1 ];
                end
                for(i_cr = 1:4)
                    kp_vals(i_cr, i, i_period) = kp_m(idx_lst(i_cr), tip_lst(i_cr), load_ind);
                    kr_vals(i_cr, i, i_period) = kr_m(idx_lst(i_cr), tip_lst(i_cr), load_ind);
                    kl_vals(i_cr, i, i_period) = kl_m(idx_lst(i_cr), tip_lst(i_cr), load_ind);

                    ep_vals(i_cr, i, i_period) = (sp_m(idx_lst(i_cr), tip_lst(i_cr), load_ind)^2 + ...
                                                  sp_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_ind)^2) / (coef(i_cr)*pi*chl);   
                    er_vals(i_cr, i, i_period) = (sr_m(idx_lst(i_cr), tip_lst(i_cr), load_ind)^2 + ...
                                                  sr_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_ind)^2) / (coef(i_cr)*pi*chl);   
                    el_vals(i_cr, i, i_period) = (sl_m(idx_lst(i_cr), tip_lst(i_cr), load_ind)^2 + ...
                                                  sl_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_ind)^2) / (coef(i_cr)*pi*chl);                
                end
            end
        end
        
        if(i_defl_type == 2)
            x_off_v_lst{5} = x_off_v_lst{5} / 2;
        end
        
        %% Plot ERR
        title_lst = [ sprintf("Fig %da. Square. ", 13-2*load_idx) ...          % rotation
                      sprintf("Fig %db. Square. ", 13-2*load_idx) ...          % offset
                      sprintf("Fig %da. Rect. ", 14-2*load_idx) ...
                      sprintf("Fig %db. Rect. ", 14-2*load_idx)];
        fn_lst = [ sprintf("plots/periods/Fig_%da", 13-2*load_idx) sprintf("plots/periods/Fig_%dc", 13-2*load_idx) ...
                   sprintf("plots/periods/Fig_%da", 14-2*load_idx) sprintf("plots/periods/Fig_%db", 14-2*load_idx) ];
        
        fig_rct = figure('Name', strcat(title_lst(i_defl_type+2), loadstr_lst(load_idx)), 'Color', 'w');
        hold on;
        box on;        
        fig_sqr = figure('Name', strcat(title_lst(i_defl_type), loadstr_lst(load_idx)), 'Color', 'w');        
        hold on;
        box on;
        
        cr_to_plot = 4 - 2*(load_idx == 1);
        for(i_period = [ 1 2 3 2 4 5 ])    
            
            if(i_defl_type == 2)
                x_v = x_off_v_lst{i_period};                
            end
            x_v_len = length(x_v);
            
            for(i_cr = 1:cr_to_plot)
                vals_p = ep_vals(i_cr, 1:x_v_len, i_period);
                vals_r = er_vals(i_cr, 1:x_v_len, i_period);
                vals_l = el_vals(i_cr, 1:x_v_len, i_period);
                if(i_defl_type == 2)
                    x_v = 0 : 0.01*chl : x_v(end);
                    vals_p = spline(x_off_v_lst{i_period}, ep_vals(i_cr, 1:x_v_len, i_period), x_v);
                    vals_r = spline(x_off_v_lst{i_period}, er_vals(i_cr, 1:x_v_len, i_period), x_v);
                    vals_l = spline(x_off_v_lst{i_period}, el_vals(i_cr, 1:x_v_len, i_period), x_v);
                end

                plot(x_v, vals_p, ...
                     'LineStyle', ls_lst(i_period), ...
                     'Color', clr_lst(i_cr), ...
                     'LineWidth', thc_lst(i_period));
%                 plot(x_v, vals_r, ...
%                      'LineStyle', ls_r, ...
%                      'Color', clr_lst(i_cr));
% 
%                 if(i_defl_type == 1)
%                     idx_lebl = 1:4:x_v_len;            % every 5 degrees
%                 else
%                     idx_lebl = 1:10:length(x_v);     
%                 end
% 
%                 plot(x_v(idx_lebl), vals_l(idx_lebl), ...
%                      'LineStyle', ls_l, ...
%                      'Color', clr_lst(i_cr), ...
%                      'Marker', 'o');
            end
            
            if(i_period == 3 || i_period == 5)
                if(i_defl_type == 2)
                    if i_period == 3
                        xticks(0 : 0.2 : 2);
                    else
                        xticks([ 0, 0.1 : 0.2 : 1.5 ]);
                    end
                else
                    xticks(defl_xtick_lst{1});
                end
                xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
                xlim([0 max(x_v)]);
                ylabel("$G(\theta_*) / G^0$", 'Interpreter', 'latex');
                
                figure(fig_rct);
            end
        end
        
        Leg1 = plot(nan, nan, 'k-');
        Leg2 = plot(nan, nan, 'k--');
        Leg3 = plot(nan, nan, 'k:');
        legend([Leg1 Leg2 Leg3], 'Conf. 2', 'Conf. 4', 'Conf. 5');
        
        figure(fig_sqr);
        
        Leg1 = plot(nan, nan, 'k--');
        Leg2 = plot(nan, nan, 'k-');
        Leg3 = plot(nan, nan, 'k:');
        legend([Leg1 Leg2 Leg3], 'Conf. 1', 'Conf. 2', 'Conf. 3');
        
%         savefig(fig_sqr, fn_lst(i_defl_type));
        savefig(fig_rct, fn_lst(i_defl_type+2));
        
%         %% Plot kinking angle
%         title_lst = [ sprintf("Fig %dc. Square. ", 13-2*load_idx) ...          % rotation
%                       sprintf("Fig %dd. Square. ", 13-2*load_idx) ...          % offset
%                       sprintf("Fig %dc. Rect. ", 14-2*load_idx) ...
%                       sprintf("Fig %dd. Rect. ", 14-2*load_idx)];
%         fn_lst = [ sprintf("plots/periods/Fig_%dc", 13-2*load_idx) sprintf("plots/periods/Fig_%dc", 13-2*load_idx) ...
%                    sprintf("plots/periods/Fig_%dd", 14-2*load_idx) sprintf("plots/periods/Fig_%dd", 14-2*load_idx) ];
%         
%         fig_off = figure('Name', strcat(title_lst(i_defl_type+2), loadstr_lst(load_idx)), 'Color', 'w');
%         hold on;
%         box on;        
%         fig_rot = figure('Name', strcat(title_lst(i_defl_type), loadstr_lst(load_idx)), 'Color', 'w');        
%         hold on;
%         box on;
%         
%         cr_to_plot = 4 - 2*(load_idx == 1);
%         for(i_period = [ 1 2 3 2 4 5 ])    
%             
%             if(i_defl_type == 2)
%                 x_v = x_off_v_lst{i_period};                
%             end
%             x_v_len = length(x_v);
%             
%             for(i_cr = 1:cr_to_plot)
%                 vals_p = kp_vals(i_cr, 1:x_v_len, i_period);
%                 vals_r = kr_vals(i_cr, 1:x_v_len, i_period);
%                 vals_l = kl_vals(i_cr, 1:x_v_len, i_period);
%                 if(i_defl_type == 2)
%                     x_v = 0 : 0.01*chl : x_v(end);
%                     vals_p = spline(x_off_v_lst{i_period}, kp_vals(i_cr, 1:x_v_len, i_period), x_v);
%                     vals_r = spline(x_off_v_lst{i_period}, kr_vals(i_cr, 1:x_v_len, i_period), x_v);
%                     vals_l = spline(x_off_v_lst{i_period}, kl_vals(i_cr, 1:x_v_len, i_period), x_v);
%                 end
% 
%                 plot(x_v, vals_p, ...
%                      'LineStyle', ls_lst(i_period), ...
%                      'Color', clr_lst(i_cr), ...
%                      'LineWidth', thc_lst(i_period));
% %                 plot(x_v, vals_r, ...
% %                      'LineStyle', ls_r, ...
% %                      'Color', clr_lst(i_cr));
% % 
% %                 if(i_defl_type == 1)
% %                     idx_lebl = 1:4:x_v_len;            % every 5 degrees
% %                 else
% %                     idx_lebl = 1:10:length(x_v);     
% %                 end
% % 
% %                 plot(x_v(idx_lebl), vals_l(idx_lebl), ...
% %                      'LineStyle', ls_l, ...
% %                      'Color', clr_lst(i_cr), ...
% %                      'Marker', 'o');
%             end
%             
%             if(i_period == 3 || i_period == 5)
%                 if(i_defl_type == 2)
%                     if i_period == 3
%                         xticks(0 : 0.2 : 2);
%                     else
%                         xticks([ 0, 0.1 : 0.2 : 1.5 ]);
%                     end
%                 else
%                     xticks(defl_xtick_lst{1});
%                 end
%                 xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
%                 xlim([0 max(x_v)]);
%                 ylabel("Kinking angle $\theta_*$ [deg]", 'Interpreter', 'latex');
%                 yticks(-90:5:90);
% 
%         %         lg_obj = legend("show");
%         %         set(lg_obj, 'Interpreter', 'latex');
%         
%                 figure(fig_off);
%             end
%         end
%         
%         Leg1 = plot(nan, nan, 'k-');
%         Leg2 = plot(nan, nan, 'k--');
%         Leg3 = plot(nan, nan, 'k:');
%         legend([Leg1 Leg2 Leg3], 'Conf. 2', 'Conf. 4', 'Conf. 5');
%         
%         figure(fig_rot);
%         
%         Leg1 = plot(nan, nan, 'k--');
%         Leg2 = plot(nan, nan, 'k-');
%         Leg3 = plot(nan, nan, 'k:');
%         legend([Leg1 Leg2 Leg3], 'Conf. 1', 'Conf. 2', 'Conf. 3');
%         
%         savefig(fig_rot, fn_lst(i_defl_type));
%         savefig(fig_off, fn_lst(i_defl_type+2));

        x_v = off_v;
    end
end

%% collinear case
function kink_postpr_for_paper_mode4()
    
    %% Common variables
    idx_lst = [ 4 3 ];                                                 
    tip_lst = [ 1 2 ];
    cr_ID_lst = [ 1 3 ];
    legstr_lst = [ "Crack 1" "Crack 3" ];
    clr_lst = [ "green" "red" ];
    ls_p = "-";
    ls_r = "--";
    ls_l = "none";    
    
    chl = 1;
    a = 2.5*chl;
    N = 4;
    load_sfx = 'ubs';
    loadstr_lst = [ "Uniaxial vertical load" "Shear load" "Hydrostatic load" ];
    ang_v = 0:90;
    off_v = 0 : 0.1*chl : (a/2);
    
    defl_xlab_lst = [ "Rotation angle $\alpha$ [deg]" "Relative uplift $\Delta / l$" ];
    defl_xtick_lst = { 0:15:90, off_v };
    
    set(0,'defaultAxesFontSize',12)
    
    %% Figure 4 and 6: SIF I, II, ERR and kinking angle, a = 2.5L
    load_idx = 1;
    plot_folder = "plots/collinear/temp";
    x_v = ang_v;
    for(i_defl_type = 1:2)
        
        % 'p' - by C.-Rice and Panasuk,'r' - by C.-Rice and maxERR, 'l' - by Leblond
        sp_vals = zeros(2, length(x_v), 2);            % SIF I and II
        ep_vals = zeros(2, length(x_v));               % ERR        
        kp_vals = zeros(2, length(x_v));               % Kinking
        sr_vals = zeros(2, length(x_v), 2);
        er_vals = zeros(2, length(x_v));
        kr_vals = zeros(2, length(x_v));                        
        sl_vals = zeros(2, length(x_v), 2);
        el_vals = zeros(2, length(x_v));
        kl_vals = zeros(2, length(x_v));                          
        for(i = 1:length(x_v))

            dn_lst = kink_create_out_dirs(false, i_defl_type, 4, chl, a, N, load_sfx);
            if(i_defl_type == 1)
                sif_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", x_v(i)));
            else
                sif_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", x_v(i)));
            end
            load(sif_fn, "sif_mat");
            
            [sp_m,kp_m] = eval_kink_panasuk(sif_mat);                           % Panasuk
            [sr_m,kr_m] = eval_kink_maxerr(sif_mat);                            % Rice
            [sl_m,kl_m] = eval_kink_lebl(sif_mat);

%             sp_m(isnan(sp_m)) = 0.0;
%             sr_m(isnan(sr_m)) = 0.0;
            kp_m(isnan(kp_m)) = 0.0;
            kr_m(isnan(kr_m)) = 0.0;
            
            for(i_cr = 1:2)
                kp_vals(i_cr, i) = kp_m(idx_lst(i_cr), tip_lst(i_cr), load_idx);
                kr_vals(i_cr, i) = kr_m(idx_lst(i_cr), tip_lst(i_cr), load_idx);
                kl_vals(i_cr, i) = kl_m(idx_lst(i_cr), tip_lst(i_cr), load_idx);
                
                sp_vals(i_cr, i, 1) = sp_m(idx_lst(i_cr), tip_lst(i_cr), load_idx) / sqrt(pi*chl);
                sr_vals(i_cr, i, 1) = sr_m(idx_lst(i_cr), tip_lst(i_cr), load_idx) / sqrt(pi*chl);
                sl_vals(i_cr, i, 1) = sl_m(idx_lst(i_cr), tip_lst(i_cr), load_idx) / sqrt(pi*chl);
                sp_vals(i_cr, i, 2) = sp_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx) / sqrt(pi*chl);
                sr_vals(i_cr, i, 2) = sr_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx) / sqrt(pi*chl);
                sl_vals(i_cr, i, 2) = sl_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx) / sqrt(pi*chl);
                
                ep_vals(i_cr, i) = (sp_m(idx_lst(i_cr), tip_lst(i_cr), load_idx)^2 + ...
                                    sp_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx)^2) / (pi*chl);   
                er_vals(i_cr, i) = (sr_m(idx_lst(i_cr), tip_lst(i_cr), load_idx)^2 + ...
                                    sr_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx)^2) / (pi*chl);   
                el_vals(i_cr, i) = (sl_m(idx_lst(i_cr), tip_lst(i_cr), load_idx)^2 + ...
                                    sl_m(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx)^2) / (pi*chl);                
            end
            
        end
        
        %% %% Plot ERR
        title_lst = [ sprintf("Fig %da. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ...
                      sprintf("Fig %db. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ];
        fn_lst = [ sprintf("%s/Fig_%da", plot_folder, 7 - load_idx), sprintf("%s/Fig_%db", plot_folder, 7 - load_idx) ];
        
        figure('Name', title_lst(i_defl_type), 'Color', 'w');
        
        hold on;
        box on;
        
        cr_to_plot = 2;
        for(i_cr = 1:cr_to_plot)
            vals_p = ep_vals(i_cr, :);
            vals_r = er_vals(i_cr, :);
            vals_l = el_vals(i_cr, :);
            if(i_defl_type == 2)
                x_v = 0 : 0.01*chl : x_v(end);
                vals_p = spline(off_v, ep_vals(i_cr, :), x_v);
                vals_r = spline(off_v, er_vals(i_cr, :), x_v);
                vals_l = spline(off_v, el_vals(i_cr, :), x_v);
            end
            
            plot(x_v, vals_p, ...
                 'LineStyle', ls_p, ...
                 'Color', clr_lst(i_cr), ...
                 'LineWidth', 1.5);
            plot(x_v, vals_r, ...
                 'LineStyle', ls_r, ...
                 'Color', clr_lst(i_cr));
            
            if(i_defl_type == 1)
                idx_lebl = 1:4:91;            % every 5 degrees
            else
                idx_lebl = 1:10:121;     
            end
            
            plot(x_v(idx_lebl), vals_l(idx_lebl), ...
                 'LineStyle', ls_l, ...
                 'Color', clr_lst(i_cr), ...
                 'Marker', 'o');
             
            Xl = xlim();
            Yl = ylim();
            xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
            ybeg = (vals_p(50) - Yl(1)) / (Yl(2) - Yl(1));
            
            annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)));
            
        end
        
        xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
        xticks(defl_xtick_lst{i_defl_type});
        xlim([0 max(x_v)]);
        ylabel("$G(\theta_*) / G^0$", 'Interpreter', 'latex');
        %yticks(-90:20:90);
        
%         lg_obj = legend("show");
%         set(lg_obj, 'Interpreter', 'latex');

        savefig(fn_lst(i_defl_type));
        saveas(gcf, fn_lst(i_defl_type), 'png');
        
        %% %% Plot SIF I
        title_lst = [ sprintf("Fig %dc. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ...
                      sprintf("Fig %dd. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ];
        fn_lst = [ sprintf("%s/Fig_%dc", plot_folder, 7 - load_idx), sprintf("%s/Fig_%dd", plot_folder, 7 - load_idx) ];
        
        figure('Name', title_lst(i_defl_type), 'Color', 'w');
        hold on;
        box on;
        
        cr_to_plot = 2;
        for(i_cr = 1:cr_to_plot)
            vals_p = sp_vals(i_cr, :, 1);
            vals_r = sr_vals(i_cr, :, 1);
            vals_l = sl_vals(i_cr, :, 1);
            if(i_defl_type == 2)
                x_v = 0 : 0.01*chl : x_v(end);
                vals_p = spline(off_v, sp_vals(i_cr, :, 1), x_v);
                vals_r = spline(off_v, sr_vals(i_cr, :, 1), x_v);
                vals_l = spline(off_v, sl_vals(i_cr, :, 1), x_v);
            end
            
            plot(x_v, vals_p, ...
                 'LineStyle', ls_p, ...
                 'Color', clr_lst(i_cr), ...
                 'LineWidth', 1.5);
            plot(x_v, vals_r, ...
                 'LineStyle', ls_r, ...
                 'Color', clr_lst(i_cr));
            
            if(i_defl_type == 1)
                idx_lebl = 1:4:91;            % every 5 degrees
            else
                idx_lebl = 1:10:121;     
            end
            
            plot(x_v(idx_lebl), vals_l(idx_lebl), ...
                 'LineStyle', ls_l, ...
                 'Color', clr_lst(i_cr), ...
                 'Marker', 'o');
             
            Xl = xlim();
            Yl = ylim();
            xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
            ybeg = (vals_p(50) - Yl(1)) / (Yl(2) - Yl(1));
            
            annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        
        xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
        xticks(defl_xtick_lst{i_defl_type});
        xlim([0 max(x_v)]);
        ylabel("$\hat{K}_I / K_I^0$", 'Interpreter', 'latex');
        %yticks(-90:20:90);
        
        savefig(fn_lst(i_defl_type));
        saveas(gcf, fn_lst(i_defl_type), 'png');
        
%         lg_obj = legend("show");
%         set(lg_obj, 'Interpreter', 'latex');
        
%         %% Plot SIF II
%         title_lst = [ "Fig 6g. Square, a = 2.5. " ...
%                       "Fig 6h. Square, a = 2.5. " ];
%         figure('Name', strcat(title_lst(i_defl_type), loadstr_lst(load_idx)), 'Color', 'w');
%         hold on;
%         
%         for(i_cr = 1:4)
%             vals1 = sr_vals(i_cr, :, 2);
%             vals2 = sl_vals(i_cr, :, 2);
%             if(i_defl_type == 2)
%                 x_v = 0 : 0.01*chl : x_v(end);
%                 vals1 = spline(off_v, sr_vals(i_cr, :, 2), x_v);
%                 vals2 = spline(off_v, sl_vals(i_cr, :, 2), x_v);
%             end
%             
%             plot(x_v, vals1, ...
%                  'LineStyle', ls_lst(i_cr), ...
%                  'Color', clr_lst(i_cr), ...
%                  'DisplayName', strcat(legstr_lst(i_cr), ". Rice"));
%             plot(x_v, vals2, ...
%                  'LineStyle', ls_lst(i_cr), ...
%                  'Color', clr_lst(i_cr), ...
%                  'LineWidth', 2, ...
%                  'DisplayName', strcat(legstr_lst(i_cr), ". Leblond"));
%         end
%         
%         xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
%         xticks(defl_xtick_lst{i_defl_type});
%         xlim([0 max(x_v)]);
%         ylabel("$\hat{K}_{II} / K_{I}^0$", 'Interpreter', 'latex');
%         %yticks(-90:20:90);
%         
%         lg_obj = legend("show");
%         set(lg_obj, 'Interpreter', 'latex');
        
        %% Plot kinking angle
        title_lst = [ sprintf("Fig %de. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ...
                      sprintf("Fig %df. Square, a = 2.5. %s", 7 - load_idx, loadstr_lst(load_idx)) ];
        fn_lst = [ sprintf("%s/Fig_%de", plot_folder, 7 - load_idx), sprintf("%s/Fig_%df", plot_folder, 7 - load_idx) ];
        
        figure('Name', title_lst(i_defl_type), 'Color', 'w');
        hold on;
        box on;
        
        cr_to_plot = 2;
        for(i_cr = 1:cr_to_plot)
            vals_p = kp_vals(i_cr, :);
            vals_r = kr_vals(i_cr, :);
            vals_l = kl_vals(i_cr, :);
            if(i_defl_type == 2)
                x_v = 0 : 0.01*chl : x_v(end);
                vals_p = spline(off_v, kp_vals(i_cr, :), x_v);
                vals_r = spline(off_v, kr_vals(i_cr, :), x_v);
                vals_l = spline(off_v, kl_vals(i_cr, :), x_v);
            end
            
            plot(x_v, vals_p, ...
                 'LineStyle', ls_p, ...
                 'Color', clr_lst(i_cr), ...
                 'LineWidth', 1.5);
            plot(x_v, vals_r, ...
                 'LineStyle', ls_r, ...
                 'Color', clr_lst(i_cr));
            
            if(i_defl_type == 1)
                idx_lebl = 1:4:91;            % every 5 degrees
            else
                idx_lebl = 1:10:121;     
            end
            
            plot(x_v(idx_lebl), vals_l(idx_lebl), ...
                 'LineStyle', ls_l, ...
                 'Color', clr_lst(i_cr), ...
                 'Marker', 'o');
             
            Xl = xlim();
            Yl = ylim();
            xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
            ybeg = (vals_p(50) - Yl(1)) / (Yl(2) - Yl(1));
            
            annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        
        xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
        xticks(defl_xtick_lst{i_defl_type});
        xlim([0 max(x_v)]);
        ylabel("Kinking angle $\theta_*$ [deg]", 'Interpreter', 'latex');
        yticks(-90:5:90);
        
%         Leg1 = plot(nan, nan, 'k-');
%         Leg2 = plot(nan, nan, 'k--');
%         Leg3 = plot(nan, nan, 'LineStyle', ls_l, 'Marker', 'o', 'Color', 'black');
%         legend([Leg1 Leg2 Leg3], 'One', 'Two', 'Three' );
        
        savefig(fn_lst(i_defl_type));
        saveas(gcf, fn_lst(i_defl_type), 'png');

        x_v = off_v;
    end
    
    
    %% Figure 3 or 5: ERR, SIF I and II for parent cracks
%     load_idx = 1;
    x_v = ang_v;
    for(i_defl_type = 1:2)
        
        s1_vals = zeros(4, length(x_v));            % SIF I                     
        s2_vals = zeros(4, length(x_v));            % SIF II
        err_vals= zeros(4, length(x_v));
        for(i = 1:length(x_v))

            dn_lst = kink_create_out_dirs(false, i_defl_type, 4, chl, a, N, load_sfx);
            if(i_defl_type == 1)
                sif_fn = strcat(dn_lst(1), sprintf("/phi=%d.mat", x_v(i)));
            else
                sif_fn = strcat(dn_lst(1), sprintf("/yoff=%.3f.mat", x_v(i)));
            end
            load(sif_fn, "sif_mat");

            for(i_cr = 1:2)
                s1_vals(i_cr, i) = sif_mat(idx_lst(i_cr), tip_lst(i_cr), load_idx) / sqrt(pi*chl);
                s2_vals(i_cr, i) = sif_mat(idx_lst(i_cr), tip_lst(i_cr)+2, load_idx) / sqrt(pi*chl);
                err_vals(i_cr, i) = (s1_vals(i_cr, i)^2 + s2_vals(i_cr, i)^2);
            end
            
        end
        
        %% %% Plot SIF I
        title_lst = [ sprintf("Fig %da. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ...
                      sprintf("Fig %db. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ];
        fn_lst = [ sprintf("%s/Fig_%da", plot_folder, 6 - load_idx), sprintf("%s/Fig_%db", plot_folder, 6 - load_idx) ];
        
        figure('Name', title_lst(i_defl_type), 'Color', 'w');
        hold on;
        box on;
        
        for(i_cr = 1:2)
            vals1 = s1_vals(i_cr, :);
            if(i_defl_type == 2)
                x_v = 0 : 0.01*chl : x_v(end);
                vals1 = spline(off_v, s1_vals(i_cr, :), x_v);
            end
            
            plot(x_v, vals1, ...
                 'LineStyle', '-', ...
                 'Color', clr_lst(i_cr), ...
                 'LineWidth', 1.5);
            
            Xl = xlim();
            Yl = ylim();
            xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
            ybeg = (vals1(50) - Yl(1)) / (Yl(2) - Yl(1));
            
            annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        
        xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
        xticks(defl_xtick_lst{i_defl_type});
        xlim([0 max(x_v)]);
        ylabel("$K_I / K_I^0$", 'Interpreter', 'latex');
        %yticks(-90:20:90);
        
%         lg_obj = legend("show");
%         set(lg_obj, 'Interpreter', 'latex');

        savefig(fn_lst(i_defl_type));
        saveas(gcf, fn_lst(i_defl_type), 'png');
        
        %% Plot SIF II
        title_lst = [ sprintf("Fig %dc. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ...
                      sprintf("Fig %dd. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ];
        fn_lst = [ sprintf("%s/Fig_%dc", plot_folder, 6 - load_idx), sprintf("%s/Fig_%dd", plot_folder, 6 - load_idx) ];
        
        figure('Name', title_lst(i_defl_type), 'Color', 'w');
        hold on;
        box on;
        
        for(i_cr = 1:2)
            vals1 = s2_vals(i_cr, :);
            if(i_defl_type == 2)
                x_v = 0 : 0.01*chl : x_v(end);
                vals1 = spline(off_v, s2_vals(i_cr, :), x_v);
            end
            
            plot(x_v, vals1, ...
                 'LineStyle', '-', ...
                 'Color', clr_lst(i_cr), ...
                 'LineWidth', 1.5);
             
            Xl = xlim();
            Yl = ylim();
            xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
            ybeg = (vals1(50) - Yl(1)) / (Yl(2) - Yl(1));
            
            annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        
        xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
        xticks(defl_xtick_lst{i_defl_type});
        xlim([0 max(x_v)]);
        ylabel("$K_{II} / K_{I}^0$", 'Interpreter', 'latex');
        %yticks(-90:20:90);
        
%         lg_obj = legend("show");
%         set(lg_obj, 'Interpreter', 'latex');

        savefig(fn_lst(i_defl_type));
        saveas(gcf, fn_lst(i_defl_type), 'png');
        
        %% Plot ERR
        title_lst = [ sprintf("Fig %de. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ...
                      sprintf("Fig %df. Square, a = 2.5. %s", 6 - load_idx, loadstr_lst(load_idx)) ];
        fn_lst = [ sprintf("%s/Fig_%de", plot_folder, 6 - load_idx), sprintf("%s/Fig_%df", plot_folder, 6 - load_idx) ];
        
        figure('Name', title_lst(i_defl_type), 'Color', 'w');
        hold on;
        box on;
        
        for(i_cr = 1:2)
            vals1 = err_vals(i_cr, :);
            if(i_defl_type == 2)
                x_v = 0 : 0.01*chl : x_v(end);
                vals1 = spline(off_v, err_vals(i_cr, :), x_v);
            end
            
            plot(x_v, vals1, ...
                 'LineStyle', '-', ...
                 'Color', clr_lst(i_cr), ...
                 'LineWidth', 1.5);
            
            Xl = xlim();
            Yl = ylim();
            xbeg = (x_v(50) - Xl(1)) / (Xl(2) - Xl(1));
            ybeg = (vals1(50) - Yl(1)) / (Yl(2) - Yl(1));
            
            annotation('textbox', [ xbeg, ybeg, 0.037, 0.06 ], 'String', int2str(cr_ID_lst(i_cr)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        
        xlabel(defl_xlab_lst(i_defl_type), 'Interpreter', 'latex');
        xticks(defl_xtick_lst{i_defl_type});
        xlim([0 max(x_v)]);
        ylabel("$G / G^0$", 'Interpreter', 'latex');
        %yticks(-90:20:90);
        
%         lg_obj = legend("show");
%         set(lg_obj, 'Interpreter', 'latex');
        
        savefig(fn_lst(i_defl_type));
        saveas(gcf, fn_lst(i_defl_type), 'png');

        x_v = off_v;
    end
end


% %% Evaluates ERR in the tip of kinking 
% function err = eval_ERR_in_kink(sif, ka)
% % sif: 
% %   2-vec with SIF I and II mode
% %%%
% 
%     k = sif(1);
%     m = sif(2);
% 
%     rt = sqrt(k^2 + 8*m^2);
%     err = 128*m^4 * (k^2 + 6*m^2 + k*rt) / (16*m^2 + (k - rt)^2)^2;
%     
% %     co2 = cos(ka/2)^2;
% %     err_trig = co2 * (k^2*co2 + m^2*(4 - 3*co2) - 2*k*m*sin(ka));
%     
% %     if(abs(err - err_trig) > 1e-6)
% %         
% %         th_lst = -pi : 0.01 : pi;
% % 
% %         c11 = (3*cos(th_lst/2) + cos(3*th_lst/2))/4;
% %         c12 =-3*(sin(th_lst/2) + sin(3*th_lst/2))/4;
% %         c21 =   (sin(th_lst/2) + sin(3*th_lst/2))/4;
% %         c22 =   (cos(th_lst/2)+3*cos(3*th_lst/2))/4;
% %         
% %         s1 = c11 * k + c12 * m;
% %         s2 = c21 * k + c22 * m;
% %         err_lst = s1.^2 + s2.^2;
% % 
% % %         co2 = cos(th_lst/2).^2;
% % %         err_lst = co2 .* (k^2*co2 + m^2*(4 - 3*co2) - 2*k*m*sin(th_lst));
% % 
% %         plot(th_lst*180/pi, err_lst / pi);
% %         
% %         error("ERR not good");
% %     end
% end
% 
% %% Evaluates SIF I in the tip of kinking 
% function s1 = eval_SIF1_in_kink(sif, ka)
% % sif: 
% %   2-vec with SIF I and II mode
% % ka:
% %   Kinking angle
% %%%
% 
%     k = sif(1);
%     m = sif(2);
% 
%     c11 = (3*cos(ka/2) + cos(3*ka/2))/4;
%     c12 =-3*(sin(ka/2) + sin(3*ka/2))/4;
% 
%     s1 = c11 * k + c12 * m;
% end

function [sn_m, kk] = eval_kink_panasuk(sif_mat)
% Evaluation of kinking angle with Cotterel-Rice SIFs and Panasuk kinking
% formula

    sif_m = sif_mat;
%     sif_m(:, [1 2], :) = max(sif_mat(:, [1 2], :), 0);
    
    sif_siz = size(sif_m);
    if(length(sif_siz) < 3)
        sif_siz(3) = 1;
    end
    
    sn_m = sif_mat;
    kk = zeros(sif_siz(1), 2, sif_siz(3));
    
    for i_load = 1:sif_siz(3)
        
        kk(:, 1, i_load) = 2*atan((sif_m(:, 1, i_load) - sqrt(...
                                   sif_m(:, 1, i_load).^2 + 8*sif_m(:, 3, i_load).^2) ...
                              ) ./ sif_m(:, 3, i_load) / 4);
        kk(:, 2, i_load) = 2*atan((sif_m(:, 2, i_load) - sqrt(...
                                   sif_m(:, 2, i_load).^2 + 8*sif_m(:, 4, i_load).^2) ...
                              ) ./ sif_m(:, 4, i_load) / 4);
    end
    
    kk(isnan(kk)) = 0.0;
    
    c11 = (3*cos(kk/2) + cos(3*kk/2))/4;
    c12 =-3*(sin(kk/2) + sin(3*kk/2))/4;
    c21 =   (sin(kk/2) + sin(3*kk/2))/4;
    c22 =   (cos(kk/2)+3*cos(3*kk/2))/4;
    
    for(i_load = 1:sif_siz(3))
    for(i_tip = 1:2)
    for(i_cr = 1:sif_siz(1))
        
        k = sif_m(i_cr, i_tip, i_load);
        m = sif_m(i_cr, i_tip+2, i_load);
        
        s1 = c11(i_cr, i_tip, i_load) * k + c12(i_cr, i_tip, i_load) * m;
        s2 = c21(i_cr, i_tip, i_load) * k + c22(i_cr, i_tip, i_load) * m;
        
        sn_m(i_cr, i_tip, i_load) = s1;
        sn_m(i_cr, i_tip+2, i_load) = s2;

    end
    end
    end
    
    kk = 180*kk/pi;
end

function [sn_m, kk] = eval_kink_maxerr(sif_mat)
% Evaluation of kinking angle with Cotterel-Rice SIFs by maximizing ERR 

    sif_m = sif_mat;
    
    sif_siz = size(sif_m);
    if(length(sif_siz) < 3)
        sif_siz(3) = 1;
    end
    
    th_lst = -pi/2 : 0.01 : pi/2;

    c11 = (3*cos(th_lst/2) + cos(3*th_lst/2))/4;
    c12 =-3*(sin(th_lst/2) + sin(3*th_lst/2))/4;
    c21 =   (sin(th_lst/2) + sin(3*th_lst/2))/4;
    c22 =   (cos(th_lst/2)+3*cos(3*th_lst/2))/4;
    
    sn_m = sif_mat;
    kk = zeros(sif_siz(1), 2, sif_siz(3));
    
    for(i_load = 1:sif_siz(3))
    for(i_tip = 1:2)
    for(i_cr = 1:sif_siz(1))
        
        k = sif_m(i_cr, i_tip, i_load);
        m = sif_m(i_cr, i_tip+2, i_load);
        
        s1 = c11 * k + c12 * m;
        s2 = c21 * k + c22 * m;
        err_lst = s1.^2 + s2.^2;
        
        [~,max_idx] = max(err_lst);
        
        sn_m(i_cr, i_tip, i_load) = s1(max_idx);
        sn_m(i_cr, i_tip+2, i_load) = s2(max_idx);
        kk(i_cr, i_tip, i_load) = th_lst(max_idx);

    end
    end
    end
    
    kk = 180*kk/pi;
end





































































