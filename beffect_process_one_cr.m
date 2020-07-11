function beffect_process_one_cr()

    bst_lst = zeros(2, 2, 3);
    bst_lst(1, 1, 1) = 1;
    bst_lst(1, 2, 2) = 1;
    bst_lst(2, 1, 2) = 1;
    bst_lst(1, 1, 3) = 1;
    bst_lst(2, 2, 3) = 1;
    
    sif_mat = crack_interact([0.0, 0.0], [ 0 ], [ 1.0 ], bst_lst);
    err = transform_sif_err(sif_mat);
    dlmwrite('results/err_single_hor.txt', err(1, 2, :), 'delimiter',' ','precision', '%13.6E');
end