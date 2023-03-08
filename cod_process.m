function [sif_mat, avg_tr] = cod_process(sif_type, vert_cr, ax, ay, chlx, chly, num_split_cr, split_type, N, bst_lst) 
%
% SIFs for kink angle for square-lattice
%
%-Brief:
%   Evaluates SIF and average tractions in configurations with multiple
%   horiznotal cracks
%
%-Params:
%---sif_type:
%       1: SIF mode I and II
%       2: SIF mode III
%
%---num_split_cr:
%       Number of horizontal cracks per cell.
%
%---vert_cr:
%       1: with vertical cracks, 0: without, only horizontal remain
%
%---split_type:
%       1: horizontal cracks are split on 'num_split_cr' cracks, which
%          stay on their lines
%       2: same, but new cracks are stacked above the original. The
%          stack is centred in the cell.
%
%---chlx:
%       HALF-LENGTH of horizontal cracks.
%
%---chly:
%       HALF-LENGTH of vertical cracks.
%
%---ax:
%       Horizontal lattice period.
%
%---ay:
%       Vertical lattice period.
%
%---N:
%       Lattice thickness (number of additional layers).
%
%---bst_lst:
%       List of boundary stress tensors.
%--------------------------------------------------------------------------
    
    num_rect_lines = 4 + 4*N;
    num_col        = num_rect_lines + 1;                                            % number of columns
    num_row        = num_rect_lines + 1;                                            % number of rows
    num_cr_v       = num_rect_lines*num_col;                                        % total number of vertical cracks
    num_cr_h       = num_rect_lines*num_row*num_split_cr;                           % total number of horizontal cracks
    
    if(vert_cr == 0)
        num_cr_v = 0;
    end
    
    if(split_type == 2)
        num_cr_h = num_col*(num_split_cr-1);
    end
    
    num_cr         = num_cr_v + num_cr_h;                                           % total number of cracks
    
    
    hl_lst = chly * ones(num_cr, 1);
    hl_lst((num_cr_v+1):end) = chlx;
    phi_lst =-pi/2 * ones(num_cr, 1);                                      % minus so that required tip was below
    phi_lst((num_cr_v+1):end) = 0;
    
    cc_lst = zeros(num_cr, 2);                                                      % cracks' centers list
    
    cr_idx = 1;
    
    if(vert_cr == 1)
        xpar = 1;
        x0 = 0;
        for i_col = 1:num_col

            ypar = 1;
            y0 = -ay/2;
            for i_row = 1:(num_row-1)

                cc_lst(cr_idx, :) = [x0, y0];            
                cr_idx = cr_idx + 1;

                y0 = y0 + ypar * ay*i_row;
                ypar = -ypar;
            end

            x0 = x0 + xpar * ax*i_col;
            xpar = -xpar;
        end

        if(cr_idx - 1 ~= num_cr_v)
            abort(); 
        end
    end
    
    ypar = 1;
    y0 = 0;
    if(split_type == 1)
        spl_cntr_lst = (ax / num_split_cr/2) : (ax / num_split_cr) : (ax - ax / num_split_cr/2);
%         spl_cntr_lst = (ax / (num_split_cr+1)) : (ax / (num_split_cr+1)) : (ax - ax / (num_split_cr+1));
    else
        spl_cntr_lst = (ay / num_split_cr) : (ay / num_split_cr) : (ay - ay / num_split_cr);
    end
    for i_row = 1:num_row
        
        if(split_type == 2 && i_row == num_row-1)
            y0 = y0 + ypar * ay*i_row;
            ypar = -ypar;
            continue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        xpar = 1;
        x0_base = -ax;
        for i_col = 1:(num_col-1)
            
            for i_cr = 1:num_split_cr
                cc_lst(cr_idx + i_cr-1, :) = [x0_base + spl_cntr_lst(i_cr), y0];            
            end
            cr_idx = cr_idx + num_split_cr;
            
            x0_base = x0_base + xpar * ax*i_col;
            xpar = -xpar;
        end
        
        y0 = y0 + ypar * ay*i_row;
        ypar = -ypar;
    end
    
% geometry testing

    fig_name = sprintf('View of configuration of size %i', num_rect_lines);
    figure('Name', fig_name);
    
    plot(cc_lst(:, 1), cc_lst(:, 2), 'o');
    hold on;

    quiver(cc_lst(:, 1), cc_lst(:, 2), hl_lst.*cos(phi_lst), hl_lst.*sin(phi_lst), 0);
    quiver(cc_lst(:, 1), cc_lst(:, 2), hl_lst.*cos(phi_lst + pi), hl_lst.*sin(phi_lst + pi), 0);
    
    for(i = 1:num_cr)
        text(cc_lst(i,1), cc_lst(i,2), string(num2str(i)));
    end
    
    axis equal;
    drawnow;
    
    vert_ind_lst = [ 
                    ...                                                         %%%%%%%%%%%%%%%% vertical cracks
                    num_row - 1 ...                                             % (0, end)
                    4, ...                                                      % (0, 2)
                    2, ...                                                      % (0, 1)
                    1, ...                                                      % (0,-1)
                    3, ...                                                      % (0,-2)
                    num_row - 2, ...                                            % (0,-end)    
                    num_row + 3, ...                                            % (N_cell_h_cr + 1, 2)
                    num_row + 1, ...                                            % (N_cell_h_cr + 1, 1)
                    num_row, ...                                                % (N_cell_h_cr + 1,-1)
                    num_row + 2, ...                                            % (N_cell_h_cr + 1,-2)
                   ];
    
    ind_lst = [
                ...                                                         %%%%%%%%%%%%%%%% horizontal cracks 
                (num_col-3)*num_split_cr + 1, ...                            % (-end, 0)
                (1:(2*num_split_cr)), ...                                    % From (-N_cell_h_cr, 0) to (+N_cell_h_cr, 0) 
                 3*num_split_cr + 1, ...                                     % (N_cell_h_cr + 2, 0)
                (1:(2*num_split_cr)) + (num_col-1)*num_split_cr, ...          % From (-N_cell_h_cr, 2) to (+N_cell_h_cr, 2) 
                (num_col+2)*num_split_cr + 1, ...                            % (N_cell_h_cr + 2, 2)
                (num_col-1)*num_split_cr ...                                 % (+end, 0)            
              ];    
    
    if(vert_cr == 1)
        ind_lst = [ vert_ind_lst, ind_lst ];
        ind_lst(11:end) = ind_lst(11:end) + num_cr_v;
    end
    
%     if(sif_type == 1)
%         [sif_mat, avg_tr] = crack_interact(cc_lst, phi_lst, hl_lst, bst_lst, ind_lst);
%     else
%         [sif_mat, avg_tr] = crack_interact_mode3(cc_lst, phi_lst, hl_lst, bst_lst, ind_lst);
%     end
end
