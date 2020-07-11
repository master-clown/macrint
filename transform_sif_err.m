function err_mat = transform_sif_err(sif_mat)
%
%-Desc:
%   Transforms SIFs to energy release rate values.
%
    err_mat = sif_mat(:, [1, 2], :).^2 + sif_mat(:, [3, 4], :).^2;
end