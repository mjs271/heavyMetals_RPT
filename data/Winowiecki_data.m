% In the cell array, the rows are:
% (1) pH, (2) Sulfate, (3) Fe, (4) Pb

% g/mol for (1) Sulfate, (2) Fe +2, (3) Pb
convert = [96.06 55.845 207.2];
% for converting from mg
convert = convert * 1e3;

% Table 22-24, p. 161-163, in Winowiecki Thesis
% Summer 2001, Trial 1 and 2 are columns 1 and 2 of cell array
% Table 35-38, p. 174-177, in Winowiecki Thesis
% Fall 2001, Trial 1 and 2 are columns 5 and 6 of cell array
% Note: pH has only one data set for each of Summer/Fall
wData = cell(4, 8);

wData{1, 1} = [7.45 7.19 7.32 7.1 7.07 7.12 7.11 7.06 7.05 7.02 7.01 7 7.03 7.01 7.02];

wData{2, 1} = [9.048 3.620 0.668 .0716 0.770 0.594 0.803 0.624 0.497 0.559 0.581 0.573 0.486 0.516 0.568 0.459];
wData{2, 2} = [1.539 1.160 0.737 0.462 0.485 0.642 0.521 0.685 0.432 0.444 0.637 0.469 0.552 0.552 0.487 0.433 0.437];
wData{3, 1} = [10.8 28.05 36.83 44.92 46.86 49.53 52.13 52.36 55.09 52.08 51.21 55.42 59.27 58.88 67.63 71.26];
wData{3, 2} = [17.77 29.32 32.38 36.12 35.96 38.12 40.44 39.97 42.7 40.6 45.71 49.98 52.74 53.69 61.92 64.67 69.84];
wData{4, 1} = [0.0196 0.1046 0.0705 0.2035 0.2059 0.1744 0.1839 0.1711 0.2234 0.1717 0.1303 0.2212 0.3086 0.2412 0.3061 0.2935];
wData{4, 2} = [0.0247 0.0875 0.0967 0.01537 0.0941 0.0764 0.0863 0.075 0.0812 0.1126 0.0571 0.1365 0.2952 0.1354 0.4166 0.2419 0.1711];

% these are the x-coordinates for the single pH values
wData{1, 3} = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 14 17 20 23 26 29] / 1e2;

% these are the x-coordinates of the two trials
wData{2, 3} = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 14 15.5 18.5 21.5 24.5 27.5 30.5] / 1e2;
wData{3, 3} = wData{2, 3};
wData{4, 3} = wData{2, 3};

wData{2, 4} = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 14 15.5 17 20 23 26 29 32] / 1e2;
wData{3, 4} = wData{2, 4};
wData{4, 4} = wData{2, 4};

wData{2, 1} = wData{2, 1} / convert(1);
wData{2, 2} = wData{2, 2} / convert(1);
wData{3, 1} = wData{3, 1} / convert(2);
wData{3, 2} = wData{3, 2} / convert(2);
wData{4, 1} = wData{4, 1} / convert(3);
wData{4, 2} = wData{4, 2} / convert(3);

wData{1, 5} = [7.17 7.2 6.91 6.88 6.85 6.78 6.83 6.81 6.81 6.83 6.89 6.89 6.86 6.87 6.88 6.89 6.87];

wData{2, 5} = [26.345 24.959 22.9 22.521 18.423 8.03 2.829 1.581 0.794 0.677 0.493 1.924];
wData{2, 6} = [21.79 9.93 1.91 0.87 0.82 0.67 0.65 0.59 0.55 0.21 0.58 0.61 0.54 0.6 0.48 0.57 0.94];
wData{3, 5} = [13.8 16.11 17.91 20.75 23.15 25 26.71 27.7 28.98 30.81 40.45 33.09 35.35 40.68 42.09];
wData{3, 6} = [22.4 25.13 27.03 32.88 32.06 33.19 36.09 39.58 45.3 49.26 59.25 62.07];
wData{4, 5} = [0.0229 0.0146 0.0125 0.0156 0.0135 0.0208 0.0104 0.0063 0.0177 0.0177 0.0249 0.0114 0.0187 0.0281 0.0114];
wData{4, 6} = [0.0218 0.0198 0.0125 0.0177 0.0146 0.0198 0.0177 0.026 0.0322 0.0488 0.2649 0.0488];

% these are the x-coordinates for the single pH values
wData{1, 7} = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 14 15.5 17 20 23 26 29 32] / 1e2;

% these are the x-coordinates of the two trials
wData{2, 7} = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 17 23 29] / 1e2;
wData{3, 7} = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 14 17 20 23 26 29] / 1e2;
wData{4, 7} = wData{3, 7};

wData{2, 8} = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 14 15.5 17 20 23 26 29 32] / 1e2;
wData{3, 8} = [0.5 2 3.5 5 6.5 8 9.5 12.5 15.5 18.5 21.5 24.5] / 1e2;
wData{4, 8} = wData{3, 8};

wData{2, 5} = wData{2, 5} / convert(1);
wData{2, 6} = wData{2, 6} / convert(1);
wData{3, 5} = wData{3, 5} / convert(2);
wData{3, 6} = wData{3, 6} / convert(2);
wData{4, 5} = wData{4, 5} / convert(3);
wData{4, 6} = wData{4, 6} / convert(3);

save('Winowiecki_data.mat', 'wData');

%%

xAxCell = {'\textbf{pH}', '\textbf{Alkalinity (mol/kgw)}', '\textbf{NO$\mathbf{_3^-}$ (mol/kgw)}', '\textbf{SO$\mathbf{_4^{-2}}$ (mol/kgw)}',...
           '\textbf{Fe$\mathbf{^{+2}}$ (mol/kgw)}', '\textbf{Pb (mol/kgw)}'};
mark = {'o', 's', 'd', '+', '^', '*'};
% mark = {'+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
% these are the axes that correspond to Arora et al.
% lim_mat = [6.0 8.0; 7e-3 1e-2; 0 1e-4; 0 5e-5; 0 1.5e-3; 1e-8 1e-5];
% these axes capture the data a little better
lim_mat = [6.5 7.5; 7e-3 1e-2; 0 8e-5; 0 6e-5; 0 1.5e-3; 1e-8 2.2e-6];

figure(3)
clf
% school monitor
set(gcf, 'Position', [2000, 5, 1200, 800])
% home monitor
% set(gcf, 'Position', [1300, 0, 1200, 900])
%     [left bottom width height]

posCell = cell(6, 1);
posCell{1} = [0.05, 0.56, 0.30, 0.43];
posCell{2} = [0.37, 0.56, 0.30, 0.43];
posCell{3} = [0.69, 0.56, 0.30, 0.43];
posCell{4} = [0.05, 0.06, 0.30, 0.43];
posCell{5} = [0.37, 0.06, 0.30, 0.43];
posCell{6} = [0.69, 0.06, 0.30, 0.43];
%     [left bottom width height]

plotNo = [1 4 5 6];
dataNo = [1 2 5 6];

colors = get(gca,'colororder');

% loop over data sets
for i = 1 : 4

%     loop over compounds
    for j = 1 : 4
        subplot('Position', posCell{plotNo(j)})
        box on
        if i == 1
            hold on
        end
        if j == 1
            X = wData{1, dataNo(i) + 2};
            if i == 1 || i == 3
                scatter(wData{j, dataNo(i)}, -X, 50, mark{i}, 'LineWidth', 2.0, 'MarkerEdgeColor', colors(i, :))
            else
                continue
            end
        else
            X = wData{j, dataNo(i) + 2};
            scatter(wData{j, dataNo(i)}, -X, 50, mark{i}, 'LineWidth', 2.0)
        end
        a = get(gca, 'YTickLabel');
        set(gca,'YTickLabel', a, 'fontsize', 12)
        xlabel(xAxCell{plotNo(j)}, 'Interpreter', 'latex', 'FontSize', 18)
        xlim(lim_mat(plotNo(j), :));
        ylim([-0.4 ; 0.0]);
        if j == 1 || j == 2
            ylabel('\textbf{Depth (m)}', 'Interpreter', 'latex', 'FontSize', 18)
        else
            set(gca,'YTickLabel',[]);
        end
        if j == 4
           set(gca, 'XScale', 'log'); 
        end
        if j == 2 && i == 4
           legend('Summer 2001, Trial 1', 'Summer 2001, Trial 2', 'Fall 2001, Trial 1', 'Fall 2001, Trial 2', 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'southeast')
        end
    end
end


