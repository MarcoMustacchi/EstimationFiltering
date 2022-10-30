%% Plot properties
ax = gca;
hold on;
ax.FontUnits = 'points';
ax.FontSize = 22;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = ':';
ax.LineWidth = 1.5;
set(findobj(gca,'type','axes'),'TickLabelInterpreter','latex'); % axis number in LaTex