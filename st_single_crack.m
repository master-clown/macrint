%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SINGLE CRACK STRESS TENSOR FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-Prefix (for all stress tensor functions): st
%
%-Params:
%---x0, y0:
%       Coordinates of the center of the crack's line.
%
%---phi0:
%       Angle between the crack's line and axis X.
%
%---len:
%       Half-length of the crack.
%
%-Tests:
%   See the function 'test_st_contour()'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st_func_mode1, st_func_mode2] = st_single_crack()

    st_func_mode1 = @st_normal_tensor;
    st_func_mode2 = @st_shear_tensor;
end

function [x_loc, y_loc] = crack_lc(x, y, x0, y0, phi0)                              % 'crack local coordinates'
%
% Desc: 
%   Converts coordinates (x, y) from global CS to local CS of the crack,
%   center of which is at (x0, y0) and which rotated on phi0 degrees
%   from the global axis X.
%
    cphi0 = cos(phi0); sphi0 = sin(phi0);
    
    x_loc = cphi0*(x-x0) + sphi0*(y-y0);
    y_loc =-sphi0*(x-x0) + cphi0*(y-y0);
end

function val = st_const_a(x, y, len)                                                % 'alpha'
% Auxiliary constant function
    
    val = (x - len).^2 + y.^2;
end

function val = st_const_b(x, y, len)                                                % 'beta'
% Auxiliary constant function
    
    val = 2*(x.^2 + y.^2 - len.^2);
end

function val = st_const_g(x, y, len)                                                % 'gamma'
% Auxiliary constant function
    
    val = (x + len).^2 + y.^2;
end

function val = st_const_d(x, y, len)                                                % 'delta'
% Auxiliary constant function
    
    val = st_const_b(x, y, len) + ...
          2*sqrt(st_const_a(x, y, len).* ...
                 st_const_g(x, y, len));
end

function val = st_const_I2(x, y, len)
% Auxiliary constant function
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(st_const_d(x,y,len));
    
    val = 4*len^2 * 1 ...
                 ./ (ar + gr + dr) ...
                 ./ dr;
end

function val = st_const_I3(x, y, len)
% Auxiliary constant function
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(st_const_d(x,y,len));
    
    val = 2*len^3 * (gr - ar) ...
                 ./ (ar.*gr.*(dr).^3);
end

function val = st_const_I4(x, y, len)
% Auxiliary constant function
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(st_const_d(x,y,len));
    
    val = 2*len^2 * (gr + ar) ...
                 ./ (ar.*gr.*(dr).^3);
end

function val = st_const_I5(x, y, len)
% Auxiliary constant function
    
    d = st_const_d(x,y,len);
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(d);
        
    val =0.5*len^3 * (d.*(gr.^3 - ar.^3) + 3*ar.*gr.*(ar+gr).^2.*(gr-ar)) ...
                  ./ ((ar.*gr).^3.*(dr).^5);
end

function val = st_const_I6(x, y, len)
% Auxiliary constant function
    
    d = st_const_d(x,y,len);
    
    % roots of a, g, d
    ar = sqrt(st_const_a(x,y,len));
    gr = sqrt(st_const_g(x,y,len));
    dr = sqrt(d);
        
    val =0.5*len^2 * (d.*(gr.^3 + ar.^3) + 3*ar.*gr.*(ar+gr).^3) ...
                  ./ ((ar.*gr).^3.*(dr).^5);
end

function val = st_normal_xx(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    
    I2 = st_const_I2(x, y, len);
    I4 = st_const_I4(x, y, len);
    I6 = st_const_I6(x, y, len);
    
    val = I2 - 8*y.^2.*(I4 - y.^2.*I6);    
end

function val = st_normal_xy(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    
    I3 = st_const_I3(x, y, len);
    I4 = st_const_I4(x, y, len);
    I5 = st_const_I5(x, y, len);
    I6 = st_const_I6(x, y, len);
    
    val = 2*(-y.*I3 + x.*y.*I4 + 4*y.^3.*(I5 - x.*I6));
end

function val = st_normal_yy(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    
    I2 = st_const_I2(x, y, len);
    I4 = st_const_I4(x, y, len);
    I6 = st_const_I6(x, y, len);
    
    val = I2 + 4*y.^2.*(I4 - 2*y.^2.*I6);    
end

function st = st_normal_tensor(x_glob, y_glob, x0, y0, phi0, len)
    
    sxx = st_normal_xx(x_glob, y_glob, x0, y0, phi0, len);
    sxy = st_normal_xy(x_glob, y_glob, x0, y0, phi0, len);
    syy = st_normal_yy(x_glob, y_glob, x0, y0, phi0, len);
    
    cphi0 = cos(phi0);
    sphi0 = sin(phi0);
    rot_mat = [ cphi0, sphi0;
               -sphi0, cphi0 ];
    
    st = rot_mat' * [ sxx, sxy;
                      sxy, syy ] * rot_mat;
end

function val = st_shear_xx(x_glob, y_glob, x0, y0, phi0, len)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    
    I3 = st_const_I3(x, y, len);
    I4 = st_const_I4(x, y, len);
    I5 = st_const_I5(x, y, len);
    I6 = st_const_I6(x, y, len);
    
    val = 2*(3*y.*I3 - 3*x.*y.*I4 - 4*y.^3.*(I5 - x.*I6));
end

function val = st_shear_xy(x_glob, y_glob, x0, y0, phi0, len)
    
    val = st_normal_xx(x_glob, y_glob, x0, y0, phi0, len);
end

function val = st_shear_yy(x_glob, y_glob, x0, y0, phi0, len)
    
    val = st_normal_xy(x_glob, y_glob, x0, y0, phi0, len);
end

function st = st_shear_tensor(x_glob, y_glob, x0, y0, phi0, len)
    
    sxx = st_shear_xx(x_glob, y_glob, x0, y0, phi0, len);
    sxy = st_shear_xy(x_glob, y_glob, x0, y0, phi0, len);
    syy = st_shear_yy(x_glob, y_glob, x0, y0, phi0, len);
    
    cphi0 = cos(phi0);
    sphi0 = sin(phi0);
    rot_mat = [ cphi0, sphi0;
               -sphi0, cphi0 ];
    
    st = rot_mat' * [ sxx, sxy;
                      sxy, syy ] * rot_mat;
end