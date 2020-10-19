


clc; close all ;
keepvars = {'ii'};
clearvars('-except', keepvars{:});

addpath('classes');
addpath('utilities');

%% Supress warnings
warning('off', 'MATLAB:nearlySingularMatrix');

%% Raw Data Log Folder
logname='ped01';
%other logs are ped02,ped03,ped04
logfolder = strcat('logs/',logname,'/');


%% Node/Network configuration
configfile = 'config/nodepositions_nesl_mobile';

%% Create Network Manager
% NetworkManager(nodeconfig, datafolder, <owr_corrections>, <twr_corrections>)
nm = NetworkManager(configfile, logfolder, 'config/antennacorrections_mocap', 'config/twrcorrections_mocap', 'config/twrcorrections_mocap' );
node_ids = nm.getNodeIds();
node_names = nm.getNodeNames();

%% Rigid body IDs
nm.setRigidBodyId('ntb-mobile', 1);

%% Black-listed nodes
%nm.blacklistNode('ntb-charlie');

%% Select reference node
nm.setReferenceNode('ntb-alpha');

%% Bootstrap the node positions and clock params before filter starts
% use first second of measurements to bootstrap node states
nm.skipTime(46);
%ped01 --> skip 46
%ped02 --> skip 15
%ped03 --> skip 15
%ped04 --> skip 15
%nm.setStaticNodesToTruePosition();
nm.bootstrapNodeStates( 4.0 );
%



%% Measurement type selection
nm.enableMessageType(nm.MSGTYPE3, true);
%nm.enableMessageType(nm.MSGTYPE1, true);
%nm.enableMessageType(nm.MSGTYPE2, true);
%nm.enableMeasurementType(nm.MEASTYPE_d, true);
%nm.enableMeasurementType(nm.MEASTYPE_r, true);
nm.enableMeasurementType(nm.MEASTYPE_R, true);
%nm.setMessageMinInterval(nm.MSGTYPE3, 2.00);
nm.enableSLATS(true);


%% Set some known anchor positions if desired
%known_anchors = {'ntb-alpha', 'ntb-bravo', 'ntb-charlie', 'ntb-delta', 'ntb-hotel', 'ntb-golf', 'ntb-echo', 'ntb-foxtrot'};
known_anchors = {};
for i=1:length(known_anchors)
    aid = known_anchors{i};
    aidx = nm.getNodeIdx( aid );
    nm.nodes{ aidx }.var_p = [0.0; 0.0; 0.0];
end

%% Process Covariances
% Process and Initial noise and state vals are determined by each node object
Q = nm.getProcessVar();
P = nm.getInitialVar();

Q2 = diag([0.1 0.1 0.1]);
P2 = 5*diag([1 1 1]);
R2 = 0.005;
%% Save as movie
SAVEMOVIE = true;
if SAVEMOVIE
    vidObj = VideoWriter('video/mocap.avi');
    vidObj.FrameRate=20;
    open(vidObj);
end

%% Position Visualization
% get current true and estimated positions
pTruStatic = nm.getTrueStaticPositions();
pEstAll = nm.getTransformedPositions();

fig = cfigure(25,25); grid on; hold on; axis equal;

htruStatic = plot3(pTruStatic(:,2),pTruStatic(:,3),pTruStatic(:,4),'^k','MarkerSize',10, 'MarkerFaceColor', 'k');
hestAll = plot3(pEstAll(:,1),pEstAll(:,2),pEstAll(:,3),'r+');
h_reach = zeros(1,1);
%herrStatic = zeros(nm.getNumStaticNodes(),1);

varscale = 1.00;
for i=1:nm.getNumNodes()
    nid = node_ids(i);
    %     % if node is static, plot error line to true position
    %     if ~nm.nodes{i}.isMobile()
    %         herrStatic(i) = plot3([pTruStatic(i,2) pEstAll(i,1)],[pTruStatic(i,3) pEstAll(i,2)],[pTruStatic(i,4) pEstAll(i,3)], nm.getPlotStyle(nid), 'Color', nm.getPlotColor(nid));
    %     end
    % if node is the reference, add text saying so
    %if nm.nodes{i}.isReference()
    %    xyz = nm.nodes{i}.getTruePosition();
    %    text( xyz(1)-0.50, xyz(2) - 0.30, xyz(3), 'Reference');
    %end
    % add text for each node
    if ~nm.nodes{i}.isMobile()
        xyz = nm.nodes{i}.getTruePosition();
        text( xyz(1)-0.50, xyz(2) + 0.50, xyz(3), nm.nodes{i}.getName() );
    end
    
    % get a node's position covariance matrix
    nidx = (i-1)*5 + 1;
    Pi = P( nidx:(nidx+2), nidx:(nidx+2) ) + 0.01*[1 0 0; 0 0 0; 0 0 1];
    %hvar(i) = plotEllipse([pEstAll(i,1); pEstAll(i,2); pEstAll(i,3)], Pi, varscale);
end


center  =zeros(3,1);
width   =7.5;%16
x_zonotope=zonotope([center,width*eye(3,3)]);
x_zonotope_e1=zonotope([center,width*eye(3,3)]);
x_zonotope_e2=zonotope([center,width*eye(3,3)]);
x_zonotope_e3=zonotope([center,width*eye(3,3)]);
hvar = zeros(1,1);

h_reach = plotZono1(x_zonotope,[1,2],'r');
h_reach_e1 = plotZono1(x_zonotope_e1,[1,2],'r');
h_reach_e2 = plotZono1(x_zonotope_e2,[1,2],'r');
h_reach_e3 = plotZono1(x_zonotope_e3,[1,2],'r');


% rigid bodies in mocap
rigid_bodies = nm.dataparser.getRigidBodyIds();
hrigidbodies = zeros( length(rigid_bodies),1 );
for i=1:length(rigid_bodies)
    hrigidbodies(i) = plot3(0, 0, 0, 'sb', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2);
end

xlim([-7 7]);
ylim([0 4]);
zlim([-7 7]);
xlabel('X Position (m)', 'FontSize',14);
ylabel('Y Position (m)', 'FontSize',14);
zlabel('Z Position (m)', 'FontSize',14);
htitle = title('Network Localization (t = 0.00s)','FontSize',16);
view(180,0);
%legend(herr, nm.getNodeNames());
drawnow;

%% Replay data and run EKF
% analysis stop time
%t_stop = 200;
t_stop = 20;
% state, time, and transformed position history
s_history = [];
p_history = {};
t_history = [];
pSupremumAll_history={};
pInfimumAll_history={};

% last global time update
meas1 = nm.getNextMeasurement();
t_start = meas1.getTime();
t_last = t_start;
k = 0;


% plotting
plot_delay = 0*0.100; % sec
plot_last = t_start;

tlast_twr = 0;
period_twr = 5.00;


MOBILEID = 10;

if strcmp(logname,'ped03') || strcmp(logname,'ped04')
    MOBILEID = 8;
end

diffEnable = 1; % enable diffusions
s2 = zeros(3,1);
%s = zeros(9*5,1);
z2vector_ent =[];

while (t_last - t_start) < t_stop
    k = k + 1;
    
    % get next measurement object
    meas = nm.getNextMeasurement();
    if isempty(meas)
        break;
    end
    walltime = meas.getTime();
    z = meas.vectorize();
    R = meas.getCovariance();
    
    
    
    % delta-time uses wallclock (desktop timestamp) for now
    dt_ref = meas.getTime() - t_last;
    t_last = walltime;
    t_history = [t_history; walltime];
    
    
    
    % configure process and measurement functions
    f = @(s) nm.processFcn(s, dt_ref);
    h = @(s) nm.measurementFcn(s, meas);
    
    
    % update state estimate
    s = nm.getState();
    % to use centralized EKF
    [s, P] = ekf(f, s, P, h, z, dt_ref*Q, R);
    nm.setState(s);
    % configure process and measurement functions
    f2 = @(s2) nm.processFcn2(s2);
    
    
    s2Transform = nm.getTransformedQuadPos(nm.nodes{9}.getStatePosition');
    for ii = 1:3
        
        h2 = @(s2) nm.measurementFcn2(s2, mod(ii,3)+1);
        
        %prepare the measurements
        htemp = zeros(1,3);
        htemp(mod(ii,3)+1) = 1;
        z2= htemp*s2Transform';
        z2vector_ent = [z2vector_ent;z2];
        yl_e1{ii} = z2;
        hl_e1{ii} = htemp;
        if htemp(2) ==1
           Rl_e1{ii} = 0.78;
        else
           Rl_e1{ii} = 0.15;
        end   
        %[s2, P2] = ekf(f2, s2, P2, h2, z2, dt_ref*Q2, R2);

    end
    x_zonotope_e1 = intersectZonoStrip1(x_zonotope_e1,hl_e1,Rl_e1,yl_e1);
    
    % update state estimate
    s = nm.getState();
    % to use centralized EKF
    [s, P] = ekf(f, s, P, h, z, dt_ref*Q, R);
    nm.setState(s);
    % configure process and measurement functions
    f2 = @(s2) nm.processFcn2(s2);
    
    
     s2Transform = nm.getTransformedQuadPos(nm.nodes{9}.getStatePosition');
    
    for ii = 1 :3
        h2 = @(s2) nm.measurementFcn2(s2, mod(ii,3)+1);
        %prepare the measurements
        htemp = zeros(1,3);
        htemp(mod(ii,3)+1) = 1;
        z2= htemp*s2Transform';
        z2vector_ent = [z2vector_ent;z2];
        yl_e2{ii} = z2;
        hl_e2{ii} = htemp;
        if htemp(2) ==1
           Rl_e2{ii} = 0.78;
        else
           Rl_e2{ii} = 0.15;
        end   
        %[s2, P2] = ekf(f2, s2, P2, h2, z2, dt_ref*Q2, R2);

    end
    x_zonotope_e2 = intersectZonoStrip1(x_zonotope_e2,hl_e2,Rl_e2,yl_e2);
    
    % update state estimate
    s = nm.getState();
    % to use centralized EKF
    [s, P] = ekf(f, s, P, h, z, dt_ref*Q, R);
    nm.setState(s);
    % configure process and measurement functions
    f2 = @(s2) nm.processFcn2(s2);
    
    
    s2Transform = nm.getTransformedQuadPos(nm.nodes{9}.getStatePosition');
    for ii = 1 :2
        h2 = @(s2) nm.measurementFcn2(s2, mod(ii,3)+1);
        
        %prepare the measurements
        htemp = zeros(1,3);
        htemp(mod(ii,3)+1) = 1;
        z2= htemp*s2Transform';
        z2vector_ent = [z2vector_ent;z2];
        yl_e3{ii} = z2;
        hl_e3{ii} = htemp;
        if htemp(2) ==1
           Rl_e3{ii} = 0.78;
        else
           Rl_e3{ii} = 0.15;
        end   
        %[s2, P2] = ekf(f2, s2, P2, h2, z2, dt_ref*Q2, R2);
        
    end
    x_zonotope_e3 = intersectZonoStrip1(x_zonotope_e3,hl_e3,Rl_e3,yl_e3);
    
    zonol{1}= x_zonotope_e1;
    zonol{2}= x_zonotope_e2;
    zonol{3}= x_zonotope_e3;
%     if k>1
%         % add aggretor zonotope after the first iteration
%         zonol{3}= x_zonotope;
%     end
    x_zonotope = andAveraging1(zonol);
    x_zonotope = zonotope([x_zonotope.center,x_zonotope.generators,Q2]);
    x_zonotope = reduce(x_zonotope,'girard',20);
    
    x_zonotope_e1 = x_zonotope;
    x_zonotope_e2 = x_zonotope;
    x_zonotope_e3 = x_zonotope;
    
    pSupremumAll_history = [ pSupremumAll_history; supremum(interval( x_zonotope))];
    pInfimumAll_history = [ pInfimumAll_history; infimum(interval( x_zonotope))];
    
    updatePlotZono(h_reach,x_zonotope,[1,3],'r');
    updatePlotZono(h_reach_e1,x_zonotope_e1,[1,3],'r');
    updatePlotZono(h_reach_e2,x_zonotope_e2,[1,3],'r');
    updatePlotZono(h_reach_e3,x_zonotope_e3,[1,3],'r');
    
    % update position estimates
    pTruStatic = nm.getTrueStaticPositions();
    %pEstAll = nm.getTransformedPositions();
    
    staticnodes = pTruStatic(1:end,2:end);
    %pEstAll = [staticnodes; s2'];
    pEstAll = [staticnodes; x_zonotope.center'];
    %pEstAll = [staticnodes; s(41:43)'];
    % update plots
    set(hestAll,'xdata',pEstAll(:,1),'ydata', pEstAll(:,2),'zdata',pEstAll(:,3));
    
    
    if walltime - plot_last >= plot_delay
        plot_last = walltime;
        % update rigid bodies
        for i=1:length(rigid_bodies)
            rb = rigid_bodies(i);
            [xyz, latency] = nm.dataparser.getMocapPos(rb, walltime);
            if latency < 0.250
                set(hrigidbodies(i), 'XData', xyz(1), 'YData', xyz(2), 'ZData', xyz(3));
            end
        end
        
        % update plot title
        tstr = sprintf('Network Localization (t = %.2fs)', (t_last - t_start));
        set(htitle, 'String', tstr);
        drawnow;
    end
    
    % append state estimate & measurement to history
    s_history = [s_history; s'];
    p_history = [p_history; pEstAll];
    
    %fprintf('t = %.2f / %.2f \n', meas.getTime()-meas1.getTime(), t_stop);
    
    if SAVEMOVIE
        f = getframe(fig);
        writeVideo(vidObj,f);
    end
    
    %pause();
end

if SAVEMOVIE
    close(vidObj);
end

% save data


saveName ='cache/temp_ent';
save(saveName, 'nm', 'k', 's_history', 'p_history','t_history','t_stop','pSupremumAll_history','pInfimumAll_history');

dlmwrite('zvector_ent.csv',z2vector_ent,'delimiter', ',', 'precision', 9)


