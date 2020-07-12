function beffect_rdif_csv()

    A = gen_rdif_matlst();
    out_rdif_matlst(A);
end

function A = gen_rdif_matlst()

    num_rdif = 6;
    num_cols = 6;
    num_crck = 4;
    
    A = zeros(num_rdif, num_cols, num_crck);
    
    f_lst = [ "results/rect2/a=2.4hl/uniax.txt", ...
              "results/rect2/a=2.4hl/biax.txt", ...
              "results/rect2/a=2.4hl/shear.txt", ...
              "results/rect2/a=4hl/uniax.txt", ...
              "results/rect2/a=4hl/biax.txt" ...
              "results/rect2/a=4hl/shear.txt" ];
               
    for i_f = 1:6
        
        mat = read_file(f_lst(i_f));
        
        for i = 1:num_rdif
            for j = 1:num_crck
                
                A(i, i_f, j) = (mat(i+1, j) - mat(i, j)) / mat(i, j);
            end
        end
    end

end

function out_rdif_matlst(A)

    f_lst = [ "results/rect2/rdif/lmost.txt", ...
              "results/rect2/rdif/lcent.txt", ...
              "results/rect2/rdif/ucent.txt", ...
              "results/rect2/rdif/umost.txt" ];
               
    for i = 1:4
        fn = f_lst{i};
        
        fd = fopen(fn, "w");
        fprintf(fd, "%7s%7s%7s%7s%7s%7s\n", "2.4u", "2.4b", "2.4s", "4u", "4b", "4s");
        fprintf(fd, "%s\n", "------------------------------------------");
        fclose(fd);
        
        dlmwrite(fn, 100*A(:, :, i), "delimiter","","precision", "%7.2f", "-append");        
    end
end

function data = read_file(fname)

    f = fopen(fname, "r");
    fgetl(f);

    data = fscanf(f, "%*d%f%f%f%f", [4, 7])';
    
    fclose(f);
end