%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST FUNCTION
% Run it to test level curves of stress fields (st_ functions).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_stress_fields(st_func)
    [st_func_mode1, st_func_mode2] = st_func();

    x0 =-2.0;
    y0 =-4.0;
    ph = pi/8;
    l = 1.0;

    % `st_func_mode*()` return tensors transformed due to rotations.
    % To see properly located stress fields we have to cancel out
    % transformations
    cph = cos(ph);
    sph = sin(ph);
    rot_mat = [ cph, sph;
               -sph, cph ];
    
    mesh = linspace(-10, 10, 1000);
    mesh_size = size(mesh, 2);
    [X, Y] = meshgrid(mesh, mesh);
    
    Z = zeros(mesh_size, mesh_size, 6);
    for i = 1:mesh_size
        for j = 1:mesh_size
            st_mode1 = rot_mat * st_func_mode1(X(i, j), Y(i, j), x0, y0, ph, l) * rot_mat';
            st_mode2 = rot_mat * st_func_mode2(X(i, j), Y(i, j), x0, y0, ph, l) * rot_mat';

            Z(i, j, 1) = st_mode1(1,1);
            Z(i, j, 2) = st_mode1(1,2);
            Z(i, j, 3) = st_mode1(2,2);
            Z(i, j, 4) = st_mode2(1,1);
            Z(i, j, 5) = st_mode2(1,2);
            Z(i, j, 6) = st_mode2(2,2);
        end
    end
    
    title_lst = { '\sigma^{n}_{xx}', '\sigma^{n}_{xy}', '\sigma^{n}_{yy}', ...
                  '\sigma^{t}_{xx}', '\sigma^{t}_{xy}', '\sigma^{t}_{yy}' };
        
    for i = 1:6
        subplot(2, 3, i);
        contour(X, Y, Z(:,:,i), [ -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6 ], ...
               'ShowText','on');
        title(title_lst(i));
        xlabel('x');
        ylabel('y');
    end
end