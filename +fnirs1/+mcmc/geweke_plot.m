
function H = geweke_plot(x)
if ~isnumeric(x)
    error('geweke diagnostic designed for numeric inputs');
end
G = table2array(fnirs1.mcmc.geweke_sequence(x));
P = size(G, 2) - 1;  % number of variables

% Precompute axis limits
xrange = [ G(1, 1), G(end, 1) ];
yrange = [ 0, max([max(G(:, 2:end)), 2.1]) ];

% Divide figure layout grid
M = floor(sqrt(P));
N = ceil(P / M);
H = subplot(M, N, 1);

for j = 1:P
    subplot(M, N, j);
    plot(G(:, 1), G(:, j + 1), 'b.', 'MarkerSize', 20);
    hold on;
    box off;
    xlim(xrange);
    ylim(yrange);
    set(gca, ...
        'FontSize', 12, ...
        'TickDir', 'out', ...
        'TickLength', 1.8 * get(gca, 'TickLength') ...
        );
    plot([G(1, 1), G(end, 1)], [2, 2], ':k');
    xlabel('Offset (Samples)');
    ylabel('Geweke Criterion');
    title("Var" + string(j));
    hold off;
end
end
