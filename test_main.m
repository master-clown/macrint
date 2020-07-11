%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% TESTING CONFIGURATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% All local functions, which actually solve %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% the problems, return RELATIVE SIFs %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_main()
    %test_2collinear_vardist([0.4, 0.2, 0.1, 0.04, 0.02], [ true, true, true ]);
    %test_2collinear_vardist(0.1 : 0.1 : 2.1, [ true, true, false ]);
    %test_2stacked_vardist(0.2 : 0.1 : 2.0, [ true, true, true]);
    test_2x2_vardist(0.1 : 0.1 : 2.0, [ true, true ]);
    %test_2collinear_rotate([0.0:pi/144:5*pi/36, 6*pi/36 : pi/36 : pi/2], [ true, true, true ]);
    %test_2collinear_translate([0.0:0.01:0.5, 0.6:0.1:2.0], [ true, true, true ]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TWO COLLINEAR IDENTICAL CRACKS UNDER NORMAL TRACTION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_2collinear_vardist(cdist_lst, output_flag_lst)
%
%-Desc:
%   Method for solving a problem for two collinear cracks under 
%   normal tractions.
%
%-Params:
%---output_flag_lst:
%       List of boolean flags for output control. By indices:
%       (1): Turn on/off table output
%       (2): Turn on/off plotting of SIF mode I in inner point
%       (3): Turn on/off plotting of SIF mode I in outer point
%

    if output_flag_lst(1)
        
        fprintf('INTERACTION OF TWO COLLINEAR IDENTICAL CRACKS UNDER NORMAL TRACTION\n');
        fprintf('%10s%15s%15s\n', ...
                'Distance', ...
                'KI/K^0 inner', ...,
                'KI/K^0 outer');
    end

    eval_num = size(cdist_lst, 2);
    sif_lst = zeros(2, eval_num);                                                   % 1st row: mode I inner, 2nd: mode I outer
    
    for i = 1:eval_num
        
        cdist = cdist_lst(i);
        sif_mat = solve_2collinear_normal(cdist, 1.0);                              % for page 286 table second arg is (1 - cdist/2)/2
        sif_lst(:, i) = sif_mat([1, 3])';
        
        if output_flag_lst(1)
            fprintf('%10.3f%15.3f%15.3f\n', cdist, sif_lst(1, i), sif_lst(2, i));
        end
    end
    
    marker_lst = [ 'o', '*' ];
    for i_data = 1:2
        if output_flag_lst(i_data+1)
            plot(cdist_lst, sif_lst(i_data, :), ...
                 'Marker', marker_lst(i_data));
            hold on;
        end
    end
    
    legend('SIF Mode I, inner point', 'SIF Mode I, outer point');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sif_mat = solve_2collinear_normal(cdist, hlen)
    %
    %-Params:
    %---cdist:
    %       Distance between the cracks.
    %
    %---len:
    %       Crack half-lengths.
    %

        dist_from_zero = cdist / 2;
        phi0 = pi/4;
        cphi0 = cos(phi0);
        sphi0 = sin(phi0);
        rot_mat = [ cphi0, sphi0;
                   -sphi0, cphi0 ];
           
        cr_center_lst = [ dist_from_zero + hlen, 0.0;
                         -dist_from_zero - hlen, 0.0 ];
        cr_phi_lst = [ phi0, phi0 ];
        cr_len_lst = [ hlen, hlen ];
        
        for i_cr = 1:2
            cr_center_lst(i_cr, :) = cr_center_lst(i_cr, :)*rot_mat;
        end
        
        ext_stress_tensor = rot_mat'* [ 0.0, 0.0;
                                        0.0, 1.0 ] * rot_mat;

        sif_mat = crack_interact(cr_center_lst, cr_phi_lst, cr_len_lst, ...
                                 ext_stress_tensor) / sqrt(pi*hlen);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% TWO STACKED IDENTICAL CRACKS UNDER NORMAL TRACTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_2stacked_vardist(cdist_lst, output_flag_lst)
%
%-Desc:
%   Method for solving a problem for two stacked cracks under 
%   normal tractions.
%
%-Params:
%---output_flag_lst:
%       List of boolean flags for output control. By indices:
%       (1): Turn on/off table output
%       (2): Turn on/off plotting of SIF mode I
%       (3): Turn on/off plotting of SIF mode II
%

    if output_flag_lst(1)
        
        fprintf('INTERACTION OF TWO STACKED IDENTICAL CRACKS UNDER NORMAL TRACTION\n');
        fprintf('%10s%15s%15s\n', ...
                'Distance', ...
                'KI/K^0', ...
                'KII/K^0');
    end

    eval_num = size(cdist_lst, 2);
    sif_lst = zeros(2, eval_num);                                        
    
    for i = 1:eval_num
        
        cdist = cdist_lst(i);
        sif_mat = solve_2stacked_vardist(cdist, 1);
        sif_lst(:, i) = sif_mat([1, 5])';
        
        if output_flag_lst(1)
            fprintf('%10.3f%15.3f%15.3f\n', cdist, ...
                    sif_lst(1, i), ...
                    sif_lst(2, i));
        end
    end
    
    marker_lst = [ 'o', '*'];
    for i_data = 1:2
        if output_flag_lst(i_data+1)
            plot(cdist_lst, sif_lst(i_data, :), ...
                 'Marker', marker_lst(i_data));
            hold on;
        end
    end
    
    legend('SIF Mode I', ...
           'SIF Mode II');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sif_mat = solve_2stacked_vardist(cdist, hlen)
    %
    %-Params:
    %---cdist:
    %       Distance between the cracks.
    %
    %---len:
    %       Crack half-lengths.
    %

        phi0 = pi/4;
        cphi0 = cos(phi0);
        sphi0 = sin(phi0);
        rot_mat = [ cphi0, sphi0;
                   -sphi0, cphi0 ];
               
        dist_from_zero = cdist / 2;
        
        cr_center_lst = [ 0.0, dist_from_zero;
                          0.0,-dist_from_zero ];
        cr_phi_lst = [ phi0, phi0 ];
        cr_len_lst = [ hlen, hlen ];

        for i_cr = 1:2
            cr_center_lst(i_cr, :) = cr_center_lst(i_cr, :)*rot_mat;
        end
        
        ext_stress_tensor = rot_mat' * [ 0.0, 0.0;
                                         0.0, 1.0 ] * rot_mat;

        sif_mat = crack_interact(cr_center_lst, cr_phi_lst, cr_len_lst, ...
                                 ext_stress_tensor) / sqrt(pi*cr_len_lst(1));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% CONFIGURATION:          ---          %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%        --- ---        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%          ---          %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_2x2_vardist(cdist_lst, output_flag_lst)
%
%-Desc:
%   Method for solving a problem for two stacked cracks under 
%   normal tractions.
%
%-Params:
%---output_flag_lst:
%       List of boolean flags for output control. By indices:
%       (1): Turn on/off table output
%       (2): Turn on/off plotting of SIF mode II
%

    if output_flag_lst(1)        
        fprintf('INTERACTION OF 2x2 CRACK SYSTEM UNDER SHEAR TRACTION\n');
        fprintf('%10s%15s\n', ...
                'Distance', ...
                'KII/K^0');
    end

    eval_num = size(cdist_lst, 2);
    sif_lst = zeros(1, eval_num);                                               
    
    for i = 1:eval_num
        
        cdist = cdist_lst(i);
        sif_mat = solve_2x2(cdist, 1);
        sif_lst(1, i) = sif_mat(14);
        
        if output_flag_lst(1)
            fprintf('%10.3f%15.3f\n', cdist, ...
                    sif_lst(1, i));
        end
    end
    
    if output_flag_lst(2)
        plot(cdist_lst, sif_lst(1, :), 'Marker', 'o'); 
    end
    
    legend('SIF Mode II');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sif_mat = solve_2x2(cdist, hlen)
    %
    %-Params:
    %---cdist:
    %       Distance between the cracks.
    %
    %---len:
    %       Crack half-lengths.
    %
    %---ttype:
    %       Uniform unit traction type. 0: normal, 1: shear. 
    %

        cr_center_lst = [ 0.0,          hlen;
                         -hlen,         0.0;
                          cdist + hlen, 0.0;
                          0.0,         -hlen ];
        cr_phi_lst = [ 0, 0, 0, 0 ];
        cr_len_lst = [ hlen, hlen, hlen, hlen ];

        ext_st = [0.0, 1.0;
                  1.0, 0.0];

        sif_mat = crack_interact(cr_center_lst, cr_phi_lst, cr_len_lst, ...
                                 ext_st) / sqrt(pi*hlen);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% TWO CRACKS, THE RIGHT ONE IS BEING ROTATED %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_2collinear_rotate(angle_lst, output_flag_lst)
%
%-Desc:
%   Method for solving a problem for two collinear cracks,
%   when the right one rotates
%
%-Params:
%---output_flag_lst:
%       List of boolean flags for output control. By indices:
%       (1): Turn on/off table output
%       (2): Turn on/off plotting of SIF mode I in inner point of left crack
%       (3): Turn on/off plotting of SIF mode II in outer point of left crack
%

    if output_flag_lst(1)
        
        fprintf('INTERACTION OF TWO COLLINEAR IDENTICAL CRACKS UNDER NORMAL TRACTION\n');
        fprintf('%10s%15s%15s\n', ...
                'Angle', ...
                'KI/K^0', ...
                'KII/K^0');
    end

    eval_num = size(angle_lst, 2);
    sif_lst = zeros(2, eval_num);                                               
    
    for i = 1:eval_num
        
        angle = angle_lst(i);
        sif_mat = solve_2collinear_rotate(angle, 1);
        sif_lst(:, i) = sif_mat([3, 7])' .* [1; -1];
        
        if output_flag_lst(1)
            fprintf('%10.3f%15.3f%15.3f\n', ...
                    180*angle/pi, ...
                    sif_lst(1, i), ...
                    sif_lst(2, i));
        end
    end
    
    angle_lst = 180/pi * angle_lst;
    
    marker_lst = [ 'o', '*' ];
    for i_data = 1:2
        if output_flag_lst(i_data+1)
            plot(angle_lst, sif_lst(i_data, :), ...
                 'Marker', marker_lst(i_data));
            yticks(0.0:0.1:1.6);
            hold on;
        end
    end
    
    legend('SIF Mode I', 'SIF Mode II');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sif_mat = solve_2collinear_rotate(angle, hlen)
    %
    %-Params:
    %---angle:
    %       Rotation angle
    %
    %---len:
    %       Crack half-lengths.
    %

        phi0 = pi/4;
        cphi0 = cos(phi0);
        sphi0 = sin(phi0);
        rot_mat = [ cphi0, sphi0;
                   -sphi0, cphi0 ];
    
        cr_center_lst = [-hlen, 0.0;
                          0.1*2*hlen + hlen, 0.0 ];
        cr_phi_lst = [ phi0, phi0 + angle ];
        cr_len_lst = [ hlen, hlen ];
        
        for i_cr = 1:2
            cr_center_lst(i_cr, :) = cr_center_lst(i_cr, :)*rot_mat;
        end
        
        ext_stress_tensor = rot_mat'* [ 0.0, 0.0;
                                        0.0, 1.0 ] * rot_mat;    

        sif_mat = crack_interact(cr_center_lst, cr_phi_lst, cr_len_lst, ...
                                 ext_stress_tensor) / sqrt(pi*hlen);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TWO CRACKS, THE RIGHT ONE IS BEING VERTICALLY TRANSLATED %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_2collinear_translate(yoff_lst, output_flag_lst)
%
%-Desc:
%   Method for solving a problem for two collinear cracks under 
%   normal tractions, when the right one moves
%
%-Params:
%---output_flag_lst:
%       List of boolean flags for output control. By indices:
%       (1): Turn on/off table output
%       (2): Turn on/off plotting of SIF mode I in inner point
%       (3): Turn on/off plotting of SIF mode I in outer point
%

    if output_flag_lst(1)
        
        fprintf('INTERACTION OF TWO COLLINEAR IDENTICAL CRACKS UNDER NORMAL TRACTION\n');
        fprintf('%10s%15s%15s\n', ...
                'Y offset', ...
                'KI/K^0', ...,
                'KII/K^0');
    end

    eval_num = size(yoff_lst, 2);
    sif_lst = zeros(2, eval_num);                                               
    
    for i = 1:eval_num
        
        yoff = yoff_lst(i);
        sif_mat = solve_2collinear_translate(yoff, 1);
        sif_lst(:, i) = sif_mat([3, 7])';
        
        if output_flag_lst(1)
            fprintf('%10.3f%15.3f%15.3f\n', yoff, ...
                    sif_lst(1, i), ...
                    sif_lst(2, i));
        end
    end
    
    marker_lst = [ 'o', '*' ];
    for i_data = 1:2
        if output_flag_lst(i_data+1)
            plot(yoff_lst, sif_lst(i_data, :), ...
                 'Marker', marker_lst(i_data));
            hold on;
        end
    end
    
    legend('SIF Mode I', 'SIF Mode II');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sif_mat = solve_2collinear_translate(yoff, hlen)
    %
    %-Params:
    %---yoff:
    %       Offset of the right crack along Y direction
    %
    %---len:
    %       Crack half-lengths.
    %
    %-Returns:
    %   Relative SIF.
    %

        dist_from_zero = hlen / 10;

        cr_center_lst = [-dist_from_zero - hlen, 0.0;
                          dist_from_zero + hlen, yoff ];
        cr_phi_lst = [ 0, 0 ];
        cr_len_lst = [ hlen, hlen ];
        ext_st = [0.0, 0.0;
                  0.0, 1.0];

        sif_mat = crack_interact(cr_center_lst, cr_phi_lst, cr_len_lst, ...
                                 ext_st) / sqrt(pi*hlen);
    end
end




































