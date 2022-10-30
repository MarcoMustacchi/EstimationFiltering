function enhancedHist(cvSignal)
    ax = gca;
    ax.FontUnits = 'points';
    ax.FontSize = 22;
    ax.Title.Interpreter = 'latex';
    % ax.Title.String = ['\textbf{' strTitle '}'];
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.String = 'Parameter''s value';
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.String = 'Count';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridLineStyle = ':';
    ax.LineWidth = 1.5;
    set(findobj(gca,'type','axes'),'TickLabelInterpreter','latex'); % axis number in LaTex
    % xlim([0 length(cvSignal)]);
    set(gcf, 'Position', get(0, 'Screensize')); % gca not working
    hold on;
    histogram(cvSignal);
    hold off;
end