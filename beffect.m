function [res_hor, res_ver] = beffect(lattice_type, num_ecl, phi, bst_lst, ccd, chl)
%
%-Desc:
%   Investigates boundary effect by evaluating SIFs of central cracks and ones near
%   the boundary. One can estimate the effect by increasing 'num_ecl' param and 
%   watching SIFs of the central cracks stopping to change.
%
%-Params:
%---lattice_type:
%       Type of periodic lattice of cracks.
%       Supported values:
%           1) 1: rectangular lattice (with possible pertubation of symmetry - 
%                 see 'phi' description below).
%               _
%              |_| - a single rectangle formed by 4 cracks
%
%           2) 2: hexagonal lattice.
%
%---num_ecl -- "number of elementary cell layers":
%       Number of additional elementary cell layers which are added to basic 
%       configuration.
%
%       Example: for rectangular lattice the basic configuration consists of cracks,
%       which form 16 rectangles 4-by-4. An elementary cell here is formed by 4 
%       rectangles 2-by-2. 'num_ecl' means how many layers of elementary cells will
%       surround the basic configuration. Here, if 'num_ecl' == 1, 20 new rectangles
%       will surround the initial 16, forming by that 36 rectangles 6-by-6.
%
%---phi:
%       Perturbation of symmetry in the lattice.
%       1) If the lattice is rectangular, horizontal cracks will rotate periodically
%          on 'phi' and minus 'phi' in alternating manner.
%       2) If the lattice is hexagonal, I have no idea.
%
%---bst -- "boundary stress tensor":
%       Stress tensor of type [ sigma_xx, sigma_xy;  
%                               sigma_xy, sigma_yy; ] -- boundary conditions.
%
%---ccd -- "crack centers distance":
%       Distance between cracks' centers.
%
%---chl -- "crack half-length":
%       ...
%
%-Returns:
%---res_hor:
%       SIFs for horizontal cracks. Every row corresponds to a certain crack,
%       namely: leftmost, left central, right central and rightmost resp.
%
%---res_ver:
%       Like for res_hor, only for vertical cracks, starting with the one at the
%       bottom.
%

    if lattice_type == 1
        [res_hor, res_ver] = beffect_rect(num_ecl, phi, bst_lst, ccd, chl);
    else
        [res_hor, res_ver] = beffect_hex(num_ecl, phi, bst_lst, ccd, chl);
    end
end


function [res_hor, res_ver] = beffect_rect(num_ecl, phi, bst_lst, ccd, chl)

    num_rect_row = 4 + 4*num_ecl;
   %num_rect     = num_rect_row^2;
    num_cr_row   = num_rect_row + 1;                                                % number of cracks in a row
    num_cr_hv    = num_rect_row*num_cr_row;                                         % total number of horizontal (or vertical) cracks
    num_cr       = 2*num_cr_hv;                                                     % total number of cracks
    
    cc_lst = zeros(num_cr, 2);                                                      % cracks' centers list
    
    for i_col = 1:(num_cr_row-1)
        for i_row = 1:num_cr_row
            
            x0 = ccd*(i_col - (num_cr_row-1)/2 - 0.5);
            y0 = ccd*(i_row - (num_cr_row+1)/2);
            
            cc_lst((i_col-1)*num_cr_row + i_row            , :) = [x0, y0];         % first half is for horizontal cracks
            cc_lst((i_col-1)*num_cr_row + i_row + num_cr_hv, :) = [y0, x0];         % second -- for vertical 
        end
    end
    
    phi_lst =-pi/2 * ones(num_cr, 1);                                               % minus so that required tip was below
    phi_lst(1:2:num_cr_hv) = phi * (-1)^num_ecl;
    phi_lst(2:2:num_cr_hv) =-phi * (-1)^num_ecl;
    
    hl_lst = chl * ones(num_cr, 1);
    
% geometry testing

%     fig_name = sprintf('View of configuration of size %i', num_ecl+1);
%     figure('Name', fig_name);
% 
%     plot(cc_lst(:, 1), cc_lst(:, 2), 'o');
%     hold on;
% %     plot(cc_lst(1:2:num_cr_hv, 1), cc_lst(1:2:num_cr_hv, 2), '*');
% %     plot(cc_lst(2:2:num_cr_hv, 1), cc_lst(2:2:num_cr_hv, 2), 'o');
% %     
%     x_off = zeros(num_cr, 1);
%     x_off(1:num_cr_hv) = chl;
%     y_off = zeros(num_cr, 1);
%     y_off((num_cr_hv+1):end) = chl;
% 
%     quiver(cc_lst(:, 1) - x_off, cc_lst(:, 2) + y_off, 2*cos(phi_lst), 2*sin(phi_lst), 0);
%     drawnow;
    
    ind_lst = [ (num_cr_row+1)/2, ...                                               % leftmost
                (num_cr_row+1)/2 + ((num_cr_row-1)/2-1)*num_cr_row, ...             % left central
                (num_cr_row+1)/2 + ((num_cr_row-1)/2-0)*num_cr_row + num_cr_hv, ... % top central
                 num_cr_hv - (num_cr_row+1)/2 + 1 + num_cr_hv ];                    % topmost
    
    sif_mat = crack_interact(cc_lst, phi_lst, hl_lst, bst_lst, ind_lst);
    
    res_hor = sif_mat(1:2, :, :);
    res_ver = sif_mat(3:4, :, :);
                    
end

function beffect_hex(num_ecl, phi, bst, ccd, chl)

end
