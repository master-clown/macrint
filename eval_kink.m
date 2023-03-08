% function res = eval_kink(sif_mat)
% % Evaluation of kinking angle with Cotterel-Rice SIFs and Panasuk kinking
% % formula
% 
%     sif_m = sif_mat;
% %     sif_m(:, [1 2], :) = max(sif_mat(:, [1 2], :), 0);
%     
%     sif_siz = size(sif_m);
%     if(length(sif_siz) < 3)
%         sif_siz(3) = 1;
%     end
%     
%     res = zeros(sif_siz(1), 2, sif_siz(3));
%     
%     for i_load = 1:sif_siz(3)
%         
%         res(:, 1, i_load) = 2*atan((sif_m(:, 1, i_load) - sqrt(...
%                                sif_m(:, 1, i_load).^2 + 8*sif_m(:, 3, i_load).^2) ...
%                               ) ./ sif_m(:, 3, i_load) / 4);
%         res(:, 2, i_load) = 2*atan((sif_m(:, 2, i_load) - sqrt(...
%                                sif_m(:, 2, i_load).^2 + 8*sif_m(:, 4, i_load).^2) ...
%                               ) ./ sif_m(:, 4, i_load) / 4);
%     end
%     
%     res = 180*res/pi;
% end

function [sn_m, kk] = eval_kink(sif_mat)
% Evaluation of kinking angle with Cotterel-Rice SIFs by maximizing ERR 

    sif_m = sif_mat;
    
    sif_siz = size(sif_m);
    if(length(sif_siz) < 3)
        sif_siz(3) = 1;
    end
    
    th_lst = -pi/2 : 0.01 : pi/2;

    c11 = (3*cos(th_lst/2) + cos(3*th_lst/2))/4;
    c12 =-3*(sin(th_lst/2) + sin(3*th_lst/2))/4;
    c21 =   (sin(th_lst/2) + sin(3*th_lst/2))/4;
    c22 =   (cos(th_lst/2)+3*cos(3*th_lst/2))/4;
    
    sn_m = sif_mat;
    kk = zeros(sif_siz(1), 2, sif_siz(3));
    
    for(i_load = 1:sif_siz(3))
    for(i_tip = 1:2)
    for(i_cr = 1:sif_siz(1))
        
        k = sif_m(i_cr, i_tip, i_load);
        m = sif_m(i_cr, i_tip+2, i_load);
        
        s1 = c11 * k + c12 * m;
        s2 = c21 * k + c22 * m;
        err_lst = s1.^2 + s2.^2;
        
        [~,max_idx] = max(err_lst);
        
        sn_m(i_cr, i_tip, i_load) = s1(max_idx);
        sn_m(i_cr, i_tip+2, i_load) = s2(max_idx);
        kk(i_cr, i_tip, i_load) = th_lst(max_idx);

    end
    end
    end
    
    kk = 180*kk/pi;
end


































