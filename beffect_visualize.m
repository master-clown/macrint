function beffect_visualize()

    hold on;
    grid on;
    
    fmt_lst_4 = { 'x-k', 's-k', 'o-k', '*-k' };
    fmt_lst_24= { 'x--k', 's--k', 'o--k', '*--k' };
    xlim([4, 28]);
    xticks(4 + 4 * (0 : 1 : 6));
    xlabel('$N$','Interpreter','Latex');
    ylabel('$G/G_0$','Interpreter','Latex');
    set(gcf,'Position',[500, 500, 500,210]);    
    %set(get(gca,'YLabel'), 'rotation', 0)
    
%     set(gcf,'Name','Uniaxial','NumberTitle','off');
%     yticks(0.85 : 0.1 : 1.35);
% %     ylim([-0.1, 4.5]);
%     vis_one('results/rect2/a=4hl/uniax.txt', fmt_lst_4, false);
%     vis_one('results/rect2/a=2.4hl/uniax.txt', fmt_lst_24, false);
%     print('results/rect2/pics/uniax.eps', '-deps');
    
%     set(gcf,'Name','Biaxial','NumberTitle','off');
%     ylim([0.46, 0.68]);
%     %yticks(0.45 : 0.025 : 0.65);
%     vis_one('results/rect2/a=4hl/biax.txt', fmt_lst_4, true);
%     vis_one('results/rect2/a=2.4hl/biax.txt', fmt_lst_24, true);
%     print('results/rect2/pics/biax.eps', '-deps');
    
    set(gcf,'Name','Shear','NumberTitle','off');
    yticks(0.75 : 0.05 : 1.05);
    ylim([0.75, 1.05]);
    vis_one('results/rect2/a=4hl/shear.txt', fmt_lst_4, true);
    vis_one('results/rect2/a=2.4hl/shear.txt', fmt_lst_24, true);
    print('results/rect2/pics/shear.eps', '-deps');
end

function vis_one(fname, out_fmt_lst, is_symm)

    S = read_file(fname);
    if is_symm
        graph_ind_lst = 2:3;% 5:(-1):4;
        num_gr = 2;
    else
        graph_ind_lst = 5:(-1):4;
        num_gr = 2;
    end
    
    for i = 1:num_gr
        plot(4 + 4*S(:, 1), S(:, graph_ind_lst(i)) / pi, out_fmt_lst{2+i});
    end
end

function data = read_file(fname)

    f = fopen(fname, 'r');
    fgetl(f);

    data = fscanf(f, '%d%f%f%f%f', [5, Inf])';
    
    fclose(f);
end