function res = eval_kink(sif_mat)

    sif_siz = size(sif_mat);
    res = zeros(sif_siz(1), 2, sif_siz(3));
    
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