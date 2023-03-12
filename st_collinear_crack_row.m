%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% COLLINEAR CRACK ROW STRESS TENSOR FUNCTIONS %%%%%%%%%%%%%%%%%%%%
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
%       Length of the crack.
%
%-Tests:
%   See the function 'test_st_contour()'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st_func_mode1, st_func_mode2] = st_collinear_crack_row(period)

    st_func_mode1 = @(x_glob, y_glob, x0, y0, phi0, hl) st_normal_tensor(x_glob, y_glob, x0, y0, phi0, hl, period);
    st_func_mode2 = @(x_glob, y_glob, x0, y0, phi0, hl) st_shear_tensor (x_glob, y_glob, x0, y0, phi0, hl, period);
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

function st = st_normal_tensor(x_glob, y_glob, x0, y0, phi0, hl, period)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    [cfDer1, cfDer2] = complexFunc(x, y, hl, period);

    imCfDer1 = imag(cfDer1);
    reCfDer2 = real(cfDer2);
    imCfDer2 = imag(cfDer2);

    sxx = -imCfDer1 - y * reCfDer2; 
    syy = -imCfDer1 + y * reCfDer2;
    sxy = y * imCfDer2;
    
    cphi0 = cos(phi0);
    sphi0 = sin(phi0);
    rot_mat = [ cphi0, sphi0;
               -sphi0, cphi0 ];
    
    st = rot_mat' * [ sxx, sxy;
                      sxy, syy ] * rot_mat;
end

function st = st_shear_tensor(x_glob, y_glob, x0, y0, phi0, hl, period)
    
    [x, y] = crack_lc(x_glob, y_glob, x0, y0, phi0);
    [cfDer1, cfDer2] = complexFunc(x, y, hl, period);

    reCfDer1 = real(cfDer1);
    imCfDer1 = imag(cfDer1);
    reCfDer2 = real(cfDer2);
    imCfDer2 = imag(cfDer2);

    sxx = 2*reCfDer1 - y * imCfDer2; 
    syy = y * imCfDer2;
    sxy = -imCfDer1 - y * reCfDer2;
    
    cphi0 = cos(phi0);
    sphi0 = sin(phi0);
    rot_mat = [ cphi0, sphi0;
               -sphi0, cphi0 ];
    
    st = rot_mat' * [ sxx, sxy;
                      sxy, syy ] * rot_mat;
end

function [cfDer1, cfDer2] = complexFunc(x, y, hl, period)

    z = x + 1i*y;
    sinZ = sin(pi/period*z);
    cosZ = cos(pi/period*z);
    sinHl = sin(pi/period*hl);

    sqrtS = sqrt((sinZ - sinHl)*(sinZ + sinHl));

    sinZbySqrtS = sinZ / sqrtS;

    cfDer1 = 1i * (1 - sinZbySqrtS);
    cfDer2 = 1i*pi/period * cosZ/sqrtS * (sinZbySqrtS^2 - 1);

end

% Below go auxiliary functions, which are not called because
% we try to accelerate calculations, so they're here for illustrative purposes 
function val = auxFuncS(z, hl, period)

    sinZ = sin(pi/period*z);
    sinHl = sin(pi/period*hl);

    val = (sinZ - sinHl)*(sinZ + sinHl);
end

function val = auxFuncS_Der1(z, hl, period)

    val = pi/period*sin(2*pi/period*z);
end

function val = complexFuncDer1(x, y, hl, period)

    z = x + 1i*y;
    sinZ = sin(pi/period*z);
    sqrtS = sqrt(auxFuncS(z, hl, period));
    
    val = -1i * (sinZ/sqrtS - 1);
end

function val = complexFuncDer2(x, y, hl, period)

    z = x + 1i*y;
    sinZ = sin(pi/period*z);
    cosZ = cos(pi/period*z);
    sqrtS = sqrt(auxFuncS(z, hl, period));
    Sder = auxFuncS_Der1(z, hl, period);
    
    val = -1i * (pi/period * cosZ/sqrtS - 1/2 * sinZ * Sder / sqrtS^3);
end









































