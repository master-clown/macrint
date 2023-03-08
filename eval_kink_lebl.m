function [sn_m, kk] = eval_kink_lebl(sif_mat)
%-Evaluation of kinking angle with Leblond formuli
%
%-Returns:
%---sn_m: matrix of new SIFs in the kinking
%---kk:   matrix of kinking angles

    persistent m F11 F12 F21 F22

    if(isempty(m))
        m(1,:) = -0.5 : (1/360) : 0.5;         % from -90 to 90 degrees
        for i = 2:20
            m(i,:) = m(1,:).^i;
        end

        F11 = 1 -     3*pi^2/8*m(2,:)  +           (pi^2-5*pi^4/128)*m(4,:)  +              (pi^2/9-11*pi^4/72+119*pi^6/15360)*m(6,:)  ...
                +      5.07790*m(8,:)  -                     2.88312*m(10,:) -                                          0.0925*m(12,:) ...
                +        2.996*m(14,:) -                       4.059*m(16,:) +                                            1.63*m(18,:) +  4.1*m(20,:);

        F12 = 0 -       3*pi/2*m(1,:)  +           (10*pi/3+pi^3/16)*m(3,:)  +               (-2*pi-133*pi^3/180+59*pi^5/1280)*m(5,:)  ...
                +    12.313906*m(7,:)  -                     7.32433*m(9,:)  +                                          1.5793*m(11,:) ...
                +       4.0216*m(13,:) -                       6.915*m(15,:) +                                            4.21*m(17,:) +  4.56*m(19,:);

        F21 = 0 +         pi/2*m(1,:)  -            (4*pi/3+pi^3/48)*m(3,:)  +               (-2*pi/3+13*pi^3/30-59*pi^5/3840)*m(5,:)  ...
                -     6.176023*m(7,:)  +                     4.44112*m(9,:)  -                                          1.5340*m(11,:) ...
                -       2.0700*m(13,:) +                       4.684*m(15,:) -                                            3.95*m(17,:) - 1.32*m(19,:);

        F22 = 1 - (4+3*pi^2/8)*m(2,:)  + (8/3+29*pi^2/18-5*pi^4/128)*m(4,:)  + (-32/15-4*pi^2/9-1159*pi^4/7200+119*pi^6/15360)*m(6,:)  ...
                +     10.58254*m(8,:)  -                     4.78511*m(10,:) -                                          1.8804*m(12,:) ...
                +        7.280*m(14,:) -                       7.591*m(16,:) +                                            0.25*m(18,:) + 12.5*m(20,:);
    end

    dim_m = size(sif_mat);
    sn_m = sif_mat;
    kk = zeros(dim_m(1), 2, dim_m(3));
    
    sif_m = sif_mat;
    
    for i_load = 1:dim_m(3)
    for i_cr = 1:dim_m(1)
        
        sif_ind = [ 1, 3 ];
        for i_tip = 1:2
           
            sif1 = F11*sif_m(i_cr, sif_ind(1), i_load) + F12*sif_m(i_cr, sif_ind(2), i_load);
            sif2 = F21*sif_m(i_cr, sif_ind(1), i_load) + F22*sif_m(i_cr, sif_ind(2), i_load);
            
            err = sif1.^2 + sif2.^2;
            [~, max_ind] = max(err);
            
            sn_m(i_cr, sif_ind(1), i_load) = sif1(max_ind);
            sn_m(i_cr, sif_ind(2), i_load) = sif2(max_ind);
            kk(i_cr, i_tip, i_load) = 180*m(1, max_ind);
            
            sif_ind = sif_ind + 1;
        end        
    end
    end
end
