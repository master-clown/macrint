function beffect_main()
    
    bst_lst = zeros(2, 2, 3);
    bst_lst(1, 1, 1) = 1;
    bst_lst(1, 2, 2) = 1;
    bst_lst(2, 1, 2) = 1;
    bst_lst(1, 1, 3) = 1;
    bst_lst(2, 2, 3) = 1;

    
    out_file_lst24 = string({ ...
                              'results/rect2/a=2.4hl/uniax.txt', ...
                              'results/rect2/a=2.4hl/shear.txt', ...
                              'results/rect2/a=2.4hl/biax.txt' ...
                          });
    out_file_lst4 = string({ ...
                              'results/rect2/a=4hl/uniax.txt', ...
                              'results/rect2/a=4hl/shear.txt', ...
                              'results/rect2/a=4hl/biax.txt' ...
                          });
    for k = 6:12
        beffect_process(k, k, 1, 2.4, 0*pi/180, bst_lst, out_file_lst24);
        beffect_process(k, k, 1, 4, 0*pi/180, bst_lst, out_file_lst4);
    end
end




















