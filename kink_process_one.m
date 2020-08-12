function [sif_mat, avg_tr] = kink_process_one(sif_type, mode, phi, yoff, chl, a, N, bst_lst)
%
% SIFs for kink angle for square-lattice
%
%-Brief:
%   Finds kinking angle of the rotated on 'phi' crack, which is located in
%   the center of a lattice of cracks.
%
%-Params:
%---sif_type:
%       1: SIF mode I and II
%       2: SIF mode III
%
%---mode:
%       1: square lattice
%       2: rectangular with period proportional to l and horL = 2 * verL.
%       3: rectangular with constant period and horL = 2 * verL.
%
%---phi:
%       Angle of the crack rotation (radians). Vectorized is possible
%
%---yoff:
%       Y-axis offset of the target crack center.
%
%---chl:
%       Crack HALF-LENGTH.
%
%---a:
%       Lattice period.
%
%---N:
%       Lattice thickness (number of additional layers).
%
%---bst_lst:
%       List of boundary stress tensors.
%--------------------------------------------------------------------------

    num_rect_row = 4 + 4*N;
   %num_rect     = num_rect_row^2;
    num_cr_row   = num_rect_row + 1;                                                % number of cracks in a row
    num_cr_hv    = num_rect_row*num_cr_row;                                         % total number of horizontal (or vertical) cracks
    num_cr       = 2*num_cr_hv;                                                     % total number of cracks
    
    cc_lst = zeros(num_cr, 2);                                                      % cracks' centers list
    
    hl_lst = chl * ones(num_cr, 1);
    
    if mode ~= 1
        hl_lst(1:num_cr_hv) = 2*hl_lst(1:num_cr_hv);
    end
    
    for i_col = 1:(num_cr_row-1)
        for i_row = 1:num_cr_row
            
            x0 = a*(i_col - (num_cr_row-1)/2 - 0.5);
            y0 = a*(i_row - (num_cr_row+1)/2);
            
            cc_lst((i_col-1)*num_cr_row + i_row            , :) = [x0, y0].*[(mode ~= 1) + 1, (mode == 3) + 1];         % first half is for horizontal cracks
            cc_lst((i_col-1)*num_cr_row + i_row + num_cr_hv, :) = [y0, x0].*[(mode ~= 1) + 1, (mode == 3) + 1];         % second -- for vertical 
        end
    end
    
    % (0,0) -- center of lattice
    ind_lst = [ ...                                                        %%%%%%%%%%%%%%%% horizontal cracks 
                (num_cr_row+1)/2, ...                                      % (-end, 0)
                (num_cr_row+1)/2 + ((num_cr_row-1)/2-1)*num_cr_row, ...    % (-1, 0)
                (num_cr_row+1)/2 + ((num_cr_row-1)/2-0)*num_cr_row+1, ...  % (+1,-2)
                (num_cr_row+1)/2 + ((num_cr_row-1)/2+0)*num_cr_row, ...    % (+1, 0)   -- the one which is rotated
                (num_cr_row+1)/2 + ((num_cr_row-1)/2+0)*num_cr_row-1, ...  % (+1, 2)
                (num_cr_row+1)/2 + ((num_cr_row-1)/2+1)*num_cr_row, ...    % (+3, 0)
                num_cr_hv - (num_cr_row+1)/2 + 1                    ...    % (+end, 0)
                ...                                                        %%%%%%%%%%%%%%%% vertical cracks
                num_cr_hv - (num_cr_row+1)/2 + 1                    ...    % (0, end)
                (num_cr_row+1)/2 + ((num_cr_row-1)/2-0)*num_cr_row, ...    % (0, 1)
                (num_cr_row+1)/2 + ((num_cr_row-1)/2-1)*num_cr_row, ...    % (0,-1)
                (num_cr_row+1)/2 + ((num_cr_row-1)/2-0)*num_cr_row+1, ...  % (2, 1)
                (num_cr_row+1)/2 + ((num_cr_row-1)/2-1)*num_cr_row+1, ...  % (2,-1)
                (num_cr_row+1)/2                                    ...  % (0,-end)                
              ];
    ind_lst(8:end) = ind_lst(8:end) + num_cr_hv;
          
    phi_lst =-pi/2 * ones(num_cr, 1);                                      % minus so that required tip was below
    phi_lst(1:num_cr_hv) = 0;
    phi_lst(ind_lst(4)) = phi_lst(ind_lst(4)) + phi;
    
    cc_lst(ind_lst(4), 2) = cc_lst(ind_lst(4), 2) + yoff;
    
% geometry testing
% 
%     fig_name = sprintf('View of configuration of size %i', N+1);
%     figure('Name', fig_name);
% 
%     plot(cc_lst(:, 1), cc_lst(:, 2), 'o');
%     hold on;
% %     plot(cc_lst(1:2:num_cr_hv, 1), cc_lst(1:2:num_cr_hv, 2), '*');
% %     plot(cc_lst(2:2:num_cr_hv, 1), cc_lst(2:2:num_cr_hv, 2), 'o');
% %     
% %     x_off = zeros(num_cr, 1);
% %     x_off(1:num_cr_hv) = 2*chl;
% %     y_off = zeros(num_cr, 1);
% %     y_off((num_cr_hv+1):end) = chl;
% 
%     %quiver(cc_lst(:, 1) - x_off, cc_lst(:, 2) + y_off, 2*hl_lst.*cos(phi_lst), 2*hl_lst.*sin(phi_lst), 0);
%     quiver(cc_lst(:, 1), cc_lst(:, 2), hl_lst.*cos(phi_lst), hl_lst.*sin(phi_lst), 0);
%     drawnow;
        
    if(sif_type == 1)
        [sif_mat, avg_tr] = crack_interact(cc_lst, phi_lst, hl_lst, bst_lst, ind_lst);
    else
        [sif_mat, avg_tr] = crack_interact_mode3(cc_lst, phi_lst, hl_lst, bst_lst, ind_lst);
    end
end