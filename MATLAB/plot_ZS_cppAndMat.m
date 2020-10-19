clear all
close all

load('CMatFiles/cppZonoStrips.mat')
%get the number of var in workspace
variablesInCurrentWorkspace = who;
numVariables = length(variablesInCurrentWorkspace);
index =1;
numofiteration = numVariables/2-1;
for i=0:numofiteration
    c_i= eval(strcat('c_',num2str(i)));
    G_i= eval(strcat('G_',num2str(i)));
    zonol{index} = zonotope([c_i,G_i]);
    sup(index,:) = supremum(interval( zonol{index}));
    infi(index,:) = infimum(interval( zonol{index}));
    cen(index,:) = c_i;
    index = index +1;
end



load('cache/paper/matZS.mat');
node_ids = nm.getNodeIds();
node_names = nm.getNodeNames();

%% Plot Path
handles = [];
t_history = nm.getAllMeasurementTimes();
p_errors = [];
rb = 1;
rb_id = 'ntb-mobile';
rb_idx = nm.getNodeIdx(rb_id);
xyzm_all = [];
xyze_all = [];
errors_all = [];
cpp_err_norm = [];
cpp_err =[];
%comparison between c++ and matlab computations
enc_flt_err=[];
time_all = [];
downsample1 = 5;
sup_all = [];
infi_all=[];
cpp_xyz_all=[];
%numofiteration = min(length(t_history),numoffiles);
for i=1:numofiteration
    if i > length(p_history)
        break;
    end
    
    t = t_history(i);
    [xyz_mocap, lat] = nm.dataparser.getMocapPos( rb, t );
   %if lat < 0.050
    xyz_est = p_history{i}(rb_idx,:);
    
    xyzm_all = [xyzm_all; xyz_mocap];
    xyze_all = [xyze_all; xyz_est];
    
    time_all = [time_all; t];
    xyz_err = xyz_mocap - xyz_est;
    errors_all = [errors_all; norm(xyz_err)];
    p_errors = [p_errors; t xyz_err];
    sup_all =[sup_all; sup(i,:)];
    infi_all=[infi_all; infi(i,:)];
    cpp_xyz_all =  [cpp_xyz_all; cen(i,:)];
    
    cpp_err_norm = [cpp_err_norm; norm(cen(i,:)-xyz_mocap)];
    cpp_err = [cpp_err; t cen(i,:)-xyz_mocap];
    enc_flt_err = [enc_flt_err;norm(cen(i,:)-xyz_est)];
  %  end
end

% a little bit of filtering for plotting
xyze_all(:,1) = medfilt1(xyze_all(:,1), 3);
xyze_all(:,2) = medfilt1(xyze_all(:,2), 3);
xyze_all(:,3) = medfilt1(xyze_all(:,3), 3);

%% Plot
% ----- path -----
cfigure(20,11);
plot3(xyzm_all(1:downsample1:end,1), xyzm_all(1:downsample1:end,2), xyzm_all(1:downsample1:end,3), 'bo', 'DisplayName', 'True Position');
hold on;
plot3(cpp_xyz_all(1:downsample1:end,1), cpp_xyz_all(1:downsample1:end,2), cpp_xyz_all(1:downsample1:end,3), 's', 'MarkerFaceColor', [0 0.75 0], 'Color', [0 0.5 0], 'DisplayName', 'Estimated Position');

% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
%legend('y1','y3')
%hl = legend(['True Position'; 'Estimated Position']);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
grid on;
%leg_str = sprintf('Estimated Position, RMSE: %.3f', sqrt( mean( errors_all ).^2 ) );
set(gca, 'View', [-17 28]);
ax = gca;
ax.FontSize = 16;
lgd= legend();
lgd.FontSize = 12;
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',305,'VerticalAlignment','middle')



% ----- norm error vs. pos -----
%subplot(2,4,[7 8]);
figure
downsample1 = 5;
tstart = time_all(1);

%plot(time_all(1:downsample:end) -tstart, medfilt1(errors_all(1:downsample:end),2),'Marker', 's')
plot(1:downsample1:length(time_all),  medfilt1(errors_all(1:downsample1:end),2),'Marker', 's')

%xlabel('Time (sec)');
xlabel('Time step');
ylabel( 'Estimation error (m)');
%ylabel(hax(2), 'Distance from centroid (m)');
ylim([0 2]);
xlim([0 length(time_all)])
grid on;
ax = gca;
ax.FontSize = 16;
set(gcf, 'Position',  [50, 50, 800, 400]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%----- floating point error

figure
ds = 1;
tstart = time_all(1);

%[hax, hl1, hl2] = plotyy(time_all(1:ds:end) -tstart, medfilt1(errors_all(1:ds:end),2), time_all(1:ds:end) -tstart, xyz_dist(1:ds:end));

plot(1:downsample1:length(time_all),  enc_flt_err(1:downsample1:end),'Marker', 's')

xlabel('Time step');
ylabel( 'Floating point error (m)');
%ylabel(hax(2), 'Distance from centroid (m)');
xlim([0 length(time_all)])
ylim([0 2e-7]);
grid on;
ax = gca;
ax.FontSize = 16;
set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%%%----------------------------%%------------------------------------%%
%the true value and the bounds on each state value
figure();
downsample1 =10;
statenum =1 ;
tstart = time_all(1);

plot(1:downsample1:length(time_all),  xyzm_all(1:downsample1:end,statenum), '-', 'Color', [0 0 1]);

hold on
plot(1:downsample1:length(time_all),  sup_all(1:downsample1:end,statenum), '-x', 'Color', [0 0 0]);
plot(1:downsample1:length(time_all),  infi_all(1:downsample1:end,statenum), '-o', 'Color', [1 0 0]);


xlabel('Time step');
ylabel('X (m)');
grid on;
xlim([0 length(time_all)])
%xlim([0 20]);
ylim([-5 1]);
legend('true value','upper bound','lower bound')
ax = gca;
ax.FontSize = 16;
set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
downsample1 =10;
statenum =2 ;
tstart = time_all(1);

plot(1:downsample1:length(time_all),  xyzm_all(1:downsample1:end,statenum), '-', 'Color', [0 0 1]);

hold on
plot(1:downsample1:length(time_all),  sup_all(1:downsample1:end,statenum), '-x', 'Color', [0 0 0]);
plot(1:downsample1:length(time_all),  infi_all(1:downsample1:end,statenum), '-o', 'Color', [1 0 0]);
xlabel('Time step');
ylabel('Y (m)');
grid on;
%xlim([0 20]);
ylim([-1 4]);
xlim([0 length(time_all)])
legend('true value','upper bound','lower bound')
ax = gca;
ax.FontSize = 16;
set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
downsample1 =10;
statenum =3 ;
tstart = time_all(1);

plot(1:downsample1:length(time_all),  xyzm_all(1:downsample1:end,statenum), '-', 'Color', [0 0 1]);

hold on
plot(1:downsample1:length(time_all),  sup_all(1:downsample1:end,statenum), '-x', 'Color', [0 0 0]);
plot(1:downsample1:length(time_all),  infi_all(1:downsample1:end,statenum), '-o', 'Color', [1 0 0]);
xlabel('Time step');
ylabel('Z (m)');
grid on;
%xlim([0 20]);
ylim([-1.5 1]);
xlim([0 length(time_all)])
legend('true value','upper bound','lower bound')
ax = gca;
ax.FontSize = 16;
set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];