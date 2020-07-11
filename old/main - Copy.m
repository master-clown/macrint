function main()
    sif = test_2colinear_normal(0.02, pi/3);
    disp(sif);
end

function sif_mat = test_2colinear_normal(cdist, phi)
%
%-Desc:
%   Solves a problem for two collinear cracks under normal tractions.
%
%-Params:
%---cdist:
%       Distance between the cracks.
%
%---phi:
%       Crack lines' angle.
%

    dist_from_zero = cdist / 2;
    cphi = cos(phi);
    sphi = sin(phi);
    
    cr_center_lst = [ (1+dist_from_zero)/2, 0.0;
                     -(1+dist_from_zero)/2, 0.0 ];
    cr_phi_lst = [ phi, phi ];
    cr_len_lst = [ (1-dist_from_zero)/2, (1-dist_from_zero)/2 ];
    
 
    rot_mat = [ cphi, sphi;
               -sphi, cphi ];
           
    for i = 1:2
        cr_center_lst(i, :) = cr_center_lst(i, :)*rot_mat;
    end
            
    ext_stress_tensor = rot_mat' * [ 0.0, 0.0;
                                     0.0, 1.0 ] * rot_mat;
    ext_st = ext_stress_tensor([1, 2, 4]);
    
    sif_mat = crack_interact(cr_center_lst, cr_phi_lst, cr_len_lst, ...
                             ext_st) / sqrt(pi*cr_len_lst(1));
end













































