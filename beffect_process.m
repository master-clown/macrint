function beffect_process(cs_beg, cs_end, l, a, phi, bst_lst, out_file_lst)
%
%-Params:
%---cs_beg and cs_end: starting and ending configuration size
%
%---a: distance between crack centers.
%
%---bst: boundary stress tensor
%

    
    % columns: cs; ERR of leftmost, left central, uppermost, upper central resp.
    
    fprintf('%s: Starting calculations for %i configurations...\n', ...
            datestr(datetime('now')), cs_end - cs_beg + 1);
    
    for i = cs_beg:cs_end
        
        i_row = i - cs_beg + 1;
        fprintf('%s: Starting the %i-th configuration.\n', datestr(datetime('now')), i_row);
        
        [sif_hor, sif_ver] = beffect(1, i, phi, bst_lst, a, l);
        err_hor = transform_sif_err(sif_hor);
        err_ver = transform_sif_err(sif_ver);
        
        fprintf('%s: Done with %i-th configuration.\n', datestr(datetime('now')), i_row);
        
        err_lst = [ err_hor(1, 2, :), err_hor(2, 2, :), ...
                    err_ver(1, 2, :), err_ver(2, 2, :) ];
        output_err(out_file_lst, 1, i, err_lst);
    end

    fprintf('%s: All tasks have been done.\n', datestr(datetime('now')));
end

function output_err(out_file_lst, mode, conf_size, err_lst)

    if mode == 0
        mode_str = 'w';
    else
        mode_str = 'a';
    end
    
    num_file = numel(out_file_lst);
    
    for i_f = 1:num_file
        fid = fopen(out_file_lst(i_f), mode_str);

        if mode == 0
            fprintf(fid, '%9s%20s%20s%20s%20s\n', 'CONF SIZE', 'ERR LEFTMOST', ...
                                                  'ERR LEFTCENT','ERR UPPERCENT', ...
                                                  'ERR UPPERMOST');
        end

        num_row = size(err_lst, 1);
        for i = 1:num_row
            fprintf(fid, '%9i%20.8E%20.8E%20.8E%20.8E\n', ...
                    conf_size, err_lst(i, 1, i_f), err_lst(i, 2, i_f), ... 
                               err_lst(i, 3, i_f), err_lst(i, 4, i_f));
        end

        fclose(fid);
    end
end
