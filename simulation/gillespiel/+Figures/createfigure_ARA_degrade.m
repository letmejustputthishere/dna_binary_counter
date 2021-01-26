function createfigure_ARA_degrade(X1, ymatrix1, title)
%CREATEFIGURE(X1, ymatrix1)
%  X1:  stairs x
%  YMATRIX1:  stairs matrix data

%  Auto-generated by MATLAB on 14-Jan-2021 20:48:48

% Create figure
figure('Name', title, 'NumberTitle', 'off');

% Create axes
axes1 = axes('Position',[0.0695652173913043 0.110179640718563 0.912374581939799 0.868263473053892]);
hold(axes1,'on');

% Create multiple lines using matrix input to stairs
stairs1 = stairs(X1,ymatrix1,'LineWidth',2);
set(stairs1(1),'DisplayName','Ara','Color',[1 0 0]);

% Create ylabel
ylabel('molecules');

% Create xlabel
xlabel('time (s)');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 28800]);
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',24,'GridAlpha',0.5,'LineStyleOrder',{'-'},'XGrid','on','XMinorGrid','on','XMinorTick','on','XScale','log','YScale','log','YGrid','on','YMinorGrid',...
    'on');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.172697682752031 0.59761081974902 0.101003344481605 0.227803738317757],'LineWidth',1,'FontSize',24,'FontName','Helvetica Neue');

