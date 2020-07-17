function kink_postpr()

    kink_sif_graph();
end

function kink_rot_graph()

    chl = 1;
    a = 3*chl;
    N = 4;
    load_sfx = 'ubs';
    
    res_fn = sprintf("results/sl_kink/kink_rot_dep/chl=%.1f a=%.1fchl N=%d load=%s/result.mat", chl, a/chl, N, load_sfx);
    if exist(res_fn, 'file') ~= 2
        fprintf("%s: '%s' does not exist.\n", datestr(datetime('now')));
        return;
    end
    
    load(res_fn, "rot_ang_lst", "kink_ang_lst");
    %kink_ang_lst = 180*kink_ang_lst/pi;
    
    plot(rot_ang_lst, kink_ang_lst);
    legend("Uniax", "Shear", "Biax");
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

function kink_eval_test()

    function res = eval_kink(sif_mat)

        sif_siz = size(sif_mat);
        res = zeros(sif_siz(1), 2, sif_mat(3));

        for k = 1:sif_siz(3)
            res(:, 1, k) = 2*atan((sif_mat(:, 1, k) - sqrt(...
                                   sif_mat(:, 1, k).^2 + 8*sif_mat(:, 3, k).^2) ...
                                  ) ./ sif_mat(:, 3, k) / 4);
            res(:, 2, k) = 2*atan((sif_mat(:, 2, k) - sqrt(...
                                   sif_mat(:, 2, k).^2 + 8*sif_mat(:, 4, k).^2) ...
                                  ) ./ sif_mat(:, 4, k) / 4);
        end

        res = 180*res/pi;
    end


end









































