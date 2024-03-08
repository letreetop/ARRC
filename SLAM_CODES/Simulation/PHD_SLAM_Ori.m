clear all;  close all;  clc;    warning('off');
    
%% System parameters
Pd = 0.99;              % Probability of detection in the field of view
Ps = 0.95;
Montecarlo = 1;     
T = 1;                % Sampling time interval
NUMpar = 100;            % Number of particles
THRESHOLDprun = 1e-5;%8e-10;   % Prun threshold
THRESHOLDmerge = 2.5;   % Merge threshold
THRESHOLDdetect = 0.45;  % Landmark detecting threshold
INTENbir = 0.05;         % Intensity of birth GMs
NUMclutter = 5;

R = diag([0.8,0.3*pi/180]).^2;    
Qbirc1 = diag([0.01,0.01].^2);
 
%% Load data
load('Environment_SLAM.mat');
NUMrobot = size(Robot_Groundtruth,3);
for r = 1:NUMrobot
    Robot_Groundtruth(3,:,r) = pi_pi(Robot_Groundtruth(3,:,r));
end
T_obs = T_total;            % Observation time
figure(1); axis([250 1250 250 1250]); hold on;  plot(Landmark_Groundtruth(:,2)',Landmark_Groundtruth(:,3)','ro','MarkerFaceColor','red','MarkerSize',3);
PREDerrortraj = zeros(NUMrobot,T_total);
for r = 1:NUMrobot
    Robot_Deadreck = zeros(3,T_total,NUMrobot);
    Robot_Deadreck(:,1,r) = Robot_Groundtruth(:,1,r);
    for t = 2:T_total
        THETA = Robot_Deadreck(3,t-1,r);
        FF = [T*cos(THETA) 0;
              T*sin(THETA) 0;
              0            T];
        Robot_Deadreck(:,t,r) = Robot_Deadreck(:,t-1,r) + FF * Robot_Odometry(:,t,r);
        PREDerrortraj(r,t) = sqrt( (Robot_Deadreck(:,t,r)-Robot_Groundtruth(:,t,r))'*(Robot_Deadreck(:,t,r)-Robot_Groundtruth(:,t,r)) );
    end
    plot(Robot_Groundtruth(1,:,r),Robot_Groundtruth(2,:,r),'k-','LineWidth',1.5);
    plot(Robot_Deadreck(1,:,r),Robot_Deadreck(2,:,r),'k--','LineWidth',1.5);
end

INTENclutter = NUMclutter / (160*pi/180) * 1e-4;
% FoV of Robots 
% General FoV
ANGLEmaxfov = 90*ones(1,5);
ANGLEminfov = -90*ones(1,5);
RANGminfov = 0*ones(1,5);
RANGmaxfov = 150*ones(1,5);
FoVANGcen = 0.5*(ANGLEminfov+ANGLEmaxfov);
ZONEang = zeros(5,10);  ZONEran = zeros(5,10);
for r = 1:NUMrobot
    for i = 1:10
        ZONEang(r,i) = ANGLEminfov(1,r) + (ANGLEmaxfov(1,r)-ANGLEminfov(1,r))/10*i;
        ZONEran(r,i) = RANGminfov(1,r) + (RANGmaxfov(1,r)-RANGminfov(1,r))/10*i;
    end
end
PDzone = [Pd,Pd,Pd,Pd,Pd,Pd,Pd,Pd,Pd,0.5*Pd];


%% Compute number of detected landmarks
[EMPTrec,LANDMARKdetpt,LANDMARKdet,MAPnumtrue] = LandmarkDetection(NUMrobot,T_obs,Robot_Groundtruth,Landmark_Groundtruth,RANGminfov,RANGmaxfov,ANGLEminfov,ANGLEmaxfov);

%% Record variables
MAPesti = cell(T_obs,NUMrobot,Montecarlo); TRAJrec = cell(NUMrobot,Montecarlo); 
PHDmap = cell(T_obs,NUMrobot,Montecarlo); OSPAmap = zeros(T_obs,NUMrobot,Montecarlo); 
TRAJerror = zeros(T_obs,NUMrobot,Montecarlo); NUMlandmark = zeros(T_obs,NUMrobot,Montecarlo);
TIMErec = zeros(1,Montecarlo);

for tt = 1:Montecarlo
    for r = 1:NUMrobot
        TRAJrec{r,tt} = zeros(3,T_obs); TRAJrec{r,tt}(:,1) = Robot_Groundtruth(:,1,r);
        MAPesti{1,r,tt} = zeros(2,0); 
    end
    for t = 2:T_obs
        for r = 1:NUMrobot
            MAPesti{t,r,tt} = zeros(2,0);
            PHDmap{t,r,tt}.inten = zeros(1,0);  PHDmap{t,r,tt}.miu = zeros(2,0);    PHDmap{t,r,tt}.cov = zeros(2,2,0);
        end
    end
end
pause(0.1);

for tt = 1:Montecarlo
    tic;
    PARrobrec = cell(1,T_obs);
    %% Generating measurements
    load('Robot_Measurement.mat');
    
    %% Initial SLAM
    PARrob = cell(NUMrobot,NUMpar);
    POSEesti = zeros(3,NUMrobot);                               % Estimated pose, used to generate adpative landmarks and output estimation
    for r = 1:NUMrobot
        xx = Robot_Groundtruth(1,1,r); yy = Robot_Groundtruth(2,1,r); theta = Robot_Groundtruth(3,1,r); theta = pi_pi(theta);
        % Initial map
        Z = Robot_Measurement{r,1};
        STATEland = zeros(2,size(Z,2));
        for i = 1:size(Z,2)
            DIS = Z(1,i);    ANG = Z(2,i);
            STATEland(:,i) = [ DIS * cos(theta + ANG) + xx;DIS * sin(theta + ANG) + yy ];
            INTENbirc1(1,i) = INTENbir;
        end
%         figure(1);  hold on;    plot(STATEland(1,:),STATEland(2,:),'k*');
        for i = 1:NUMpar
            PARrob{r,i}.traj = Robot_Groundtruth(:,1,r);     % pose of robot
            PARrob{r,i}.wei = 1/NUMpar;        % weight of particle
            %-----------------------------------
            PARrob{r,i}.intenc1 = INTENbirc1;                     % weight of map GM C1
            PARrob{r,i}.miuc1 = STATEland;                       % mean vector of map GM C1
            PARrob{r,i}.covc1 = repmat(Qbirc1,1,1,size(STATEland,2));                     % covariance of map GM C1
        end
    end
    
    %% SLAM begins
    for t = 2:T_obs
        fprintf('Montecarlo: %d, Time: %d\n',tt,t);
        % Local SLAM
        MARKresample = zeros(1,NUMrobot);   BIRmea = cell(1,NUMrobot);  
        for r = 1:NUMrobot
            Qtraj = Qodo(:,:,r);
            ZM = Robot_Measurement{r,t};    
            PARlikeli = zeros(1,NUMpar); % Likelihood just computed based on class 1
            for i = 1:NUMpar
                %% Prediction for trajectory
                theta = PARrob{r,i}.traj(3,1);      % Orientation of last time
                ODOrec = Robot_Odometry(:,t,r);         % Odometry at this time
                F = [T*cos(theta),0;
                     T*sin(theta),0;
                     0,           T];
                PREDodo = PARrob{r,i}.traj + F*ODOrec;
                if isempty(ZM) == 1
                    PARrob{r,i}.traj = PREDodo;
                else
                    PARrob{r,i}.traj = PARrob{r,i}.traj + F * (ODOrec + sqrt(Qtraj)*randn(2,1));
                end
                xx = PARrob{r,i}.traj(1,1);     yy = PARrob{r,i}.traj(2,1);     theta = PARrob{r,i}.traj(3,1);
                %% Prediction of map
                % Step 1: Find the Gaussian of class 1 inside the FoV
                MIU = PARrob{r,i}.miuc1;  IDc1fov = zeros(1,0);  
                for j = 1:size(MIU,2)     % Find landmarks that inside the FoV and then attach with Pd
                    zx = MIU(1,j); zy = MIU(2,j);
                    ANGpred = atan2( zy-yy,zx-xx ) - theta;  ANGpred = pi_pi(ANGpred);  DISpred = sqrt( (zy-yy)^2+(zx-xx)^2 ); 
                    if DISpred <= RANGmaxfov(1,r)+0.6 && ANGpred*180/pi >= (ANGLEminfov(1,r))-3 && ANGpred*180/pi <= (ANGLEmaxfov(1,r))+3   % In the current FoV
                        IDc1fov = cat(2,IDc1fov,j);    
%                         PDtem = PDcompute(ZONEran(r,:),FoVANGcen(1,r),ZONEang(r,:),PDzone,DISpred,ANGpred);    GMPdc1(1,j) = PDtem;
                    end
                end
                IDc1out = setdiff([1:size(MIU,2)],IDc1fov);     % Static landmarks that outside the FoV, which not be involved in the current iteration
                INTENc1out = PARrob{r,i}.intenc1(1,IDc1out);    MIUc1out = PARrob{r,i}.miuc1(:,IDc1out);    COVc1out = PARrob{r,i}.covc1(:,:,IDc1out);
                % Landmarks that remain static
                INTENc1 = PARrob{r,i}.intenc1(1,IDc1fov);     MIUc1 = PARrob{r,i}.miuc1(:,IDc1fov);   COVc1 = PARrob{r,i}.covc1(:,:,IDc1fov);
                for j = 1:size(COVc1,3)
                    COVc1(:,:,j) = COVc1(:,:,j) + diag([0.01,0.01]);
                end
                
                
                %% Update
                if isempty(IDc1fov) ~= 1            % Update carry out if and only if there are static landmarks inside the current FoV
                    % Step 1: Static map selection
                    NUMgm = size(INTENc1,2);    Zpre = zeros(2,NUMgm);  GMPD = zeros(1,NUMgm);  GMid = zeros(1,0);
                    HH = zeros(2,2,NUMgm);      SS = zeros(2,2,NUMgm);  KK = zeros(2,2,NUMgm);  COVup = zeros(2,2,NUMgm);
                    for j = 1:NUMgm
                        zx = MIUc1(1,j); zy = MIUc1(2,j); 
                        ANGpred = atan2( zy-yy,zx-xx ) - theta;  ANGpred = pi_pi(ANGpred);  DISpred = sqrt( (zy-yy)^2+(zx-xx)^2 );     Zpre(:,j) = [DISpred;ANGpred];
                        rsqu = (zy-yy)^2+(zx-xx)^2;
                        HH(:,:,j) = [ (zx-xx)/sqrt(rsqu),(zy-yy)/sqrt(rsqu);
                                      -(zy-yy)/rsqu,     (zx-xx)/rsqu ];
                        SS(:,:,j) = HH(:,:,j) * COVc1(:,:,j) * HH(:,:,j)' + R;  KK(:,:,j) = COVc1(:,:,j) * HH(:,:,j)' * pinv(SS(:,:,j));
                        COVup(:,:,j) = (eye(2)-KK(:,:,j)*HH(:,:,j)) * COVc1(:,:,j);     % Used only for record
                        if DISpred <= RANGmaxfov(1,r) && ANGpred*180/pi >= (ANGLEminfov(1,r)) && ANGpred*180/pi <= (ANGLEmaxfov(1,r))   % In the current FoV
                            GMid = cat(2,GMid,j);
                            PDtem = PDcompute(ZONEran(r,:),FoVANGcen(1,r),ZONEang(r,:),PDzone,DISpred,ANGpred);    GMPD(1,j) = PDtem;
                        end
                    end
                    % Step 2A: Compute particle likelihood
                    LIKELIpz = zeros(1,size(ZM,2));
                    for j = 1:size(ZM,2)
                        if isempty(GMid) ~= 1
                            LIKELIpf = zeros(1,NUMgm);
                            for k = 1:NUMgm
                                zx = MIUc1(1,k); zy = MIUc1(2,k); 
                                zdiff = ZM(1:2,j) - Zpre(:,k); zdiff(2,1) = pi_pi(zdiff(2,1));
                                LIKELIpf(1,k) = 1/sqrt(det(2*pi*R))*exp( -0.5*(zdiff)'*pinv(R)*zdiff ) * INTENc1(1,k);
                            end
                            LIKELIpz(1,j) = INTENclutter + sum(LIKELIpf,2);
                            MARKresample(1,r) = 1;
                        else
                            LIKELIpz(1,j) = INTENclutter;
                        end
                    end
                    PARlikeli(1,i) = prod(LIKELIpz,2) + 1e-99;
                    % Step 2B: Update static map
                    LIKELIall = zeros(1,size(ZM,2));        % Used to normalize PHD
                    INTENupc1 = zeros(1,0);         MIUupc1 = zeros(2,0);   COVupc1 = zeros(2,2,0);   
                    INTENudc1 = (1-GMPD).*INTENc1;  MIUudc1 = MIUc1;        COVudc1 = COVc1;
                    for j = 1:size(ZM,2)
                        LIKELIgm = zeros(1,NUMgm);  zreal = ZM(1:2,j);
                        for k = 1:NUMgm
                            zpred = Zpre(:,k); zdiff = zreal-zpred; zdiff(2,1) = pi_pi(zdiff(2,1));
                            LIKELIgm(1,k) = GMPD(1,k) * (det(2*pi*SS(:,:,k)))^(-1) * exp( -0.5*zdiff'*pinv(SS(:,:,k))*zdiff ) * INTENc1(1,k) + 1e-299;
                            INTENupc1 = cat(2,INTENupc1,LIKELIgm(1,k));
                            COVupc1 = cat(3,COVupc1,COVup(:,:,k));
                            MIUupc1 = cat(2,MIUupc1, MIUc1(:,k) + KK(:,:,k)*zdiff );
                        end
                        LIKELIall(1,j) = sum(LIKELIgm,2);
                    end
                    LIKELIcomb = sum(LIKELIall,1);   
                    
                    % Step 3: Compute the final PHD of each class
                    % Step 3A: Class 1
                    NUMgm = size(INTENc1,2);
                    for j = 1:size(ZM,2)
                        Id = (j-1)*NUMgm+1:j*NUMgm;
                        INTENupc1(1,Id) = INTENupc1(1,Id) ./(INTENclutter+LIKELIcomb(1,j));
                    end
                    INTENupc1 = min(INTENupc1,0.995);
                    INTENupc1 = cat(2,INTENudc1,INTENupc1);   MIUupc1 = cat(2,MIUudc1,MIUupc1);   COVupc1 = cat(3,COVudc1,COVupc1);
                    [w1,x1,P1] = gaus_prune(INTENupc1,MIUupc1,COVupc1,THRESHOLDprun);  [w2,x2,P2] = gaus_merge_static(w1',x1,P1,THRESHOLDmerge);  
                    WEImerge = w2'; MIUmerge = x2;  COVmerge = P2;  WEImerge = min(WEImerge,0.999);
                    % Combine static landmarks that inside/outside the FoV together!
                    INTENupc1 = cat(2,INTENc1out,WEImerge);   MIUupc1 = cat(2,MIUc1out,MIUmerge);   COVupc1 = cat(3,COVc1out,COVmerge); 
                   
                    % Step 4: Copy everything into the particle
                    PARrob{r,i}.intenc1 = INTENupc1;    PARrob{r,i}.miuc1 = MIUupc1;    PARrob{r,i}.covc1 = COVupc1;
                else
                    PARrob{r,i}.traj = PREDodo;     % If there is no detected landmark inside the FoV, then trust the odometry
                    PARlikeli(1,i) = 1e-99;
                end
            end
            %% State estimation
            % Robot trajectory estimation
            MAXid = find(PARlikeli == max(PARlikeli));  
            POSEesti(:,r) = PARrob{r,MAXid(1,1)}.traj;         % Estimated pose (local coordinate)
            TRAJrec{r,tt}(:,t) = POSEesti(:,r);
            TRAJerror(t,r,tt) = sqrt( (TRAJrec{r,tt}(1:2,t) - Robot_Groundtruth(1:2,t,r))' * (TRAJrec{r,tt}(1:2,t) - Robot_Groundtruth(1:2,t,r)) );     % Pose estimation error
            % Landmark estimation
            IDmap = find(PARrob{r,MAXid(1,1)}.intenc1 >= THRESHOLDdetect);  
            PHDmap{t,r,tt}.inten = PARrob{r,MAXid(1,1)}.intenc1;
            PHDmap{t,r,tt}.miu = PARrob{r,MAXid(1,1)}.miuc1;
            PHDmap{t,r,tt}.cov = PARrob{r,MAXid(1,1)}.covc1;
            NUMlandmark(t,r,tt) = size(IDmap,2);    % Number of detected landmarks
            MAPdec = PARrob{r,MAXid(1,1)}.miuc1(:,IDmap);    
            MAPesti{t,r,tt} = MAPdec;
            MAPdetetrue = [Landmark_Groundtruth(LANDMARKdet{r,t},2)';Landmark_Groundtruth(LANDMARKdet{r,t},3)'];
            OSPAmap(t,r,tt) = ospa_dist(MAPdetetrue,MAPesti{t,r,tt},20,2);        % Map estimation error
%             figure(1);  hold on;  plot(MAPesti{t,r,tt}(1,:),MAPesti{t,r,tt}(2,:),'b*');
            
            %% Birth PHD generate
            ZM = Robot_Measurement{r,t};
            INTENbirc1 = zeros(1,0);    MIUbirc1 = zeros(2,0);  COVbirc1 = zeros(2,2,0);
            for j = 1:size(ZM,2)
                xx = POSEesti(1,r);     yy = POSEesti(2,r);     theta = POSEesti(3,r);
                DIS = ZM(1,j);    ANG = ZM(2,j);
                zx = DIS * cos(theta + ANG) + xx;   zy = DIS * sin(theta + ANG) + yy;
                MIUtot = PARrob{r,MAXid(1,1)}.miuc1;
                MATRIXdis = repmat([zx;zy],1,size(MIUtot,2)) - MIUtot;
                DISall = sqrt( MATRIXdis(1,:).^2+MATRIXdis(2,:).^2 );
                if min(DISall) >= 5     % If the measurement is not close to any existing landmark/target
                    INTENbirc1 = cat(2,INTENbirc1,INTENbir);
                    MIUbirc1 = cat(2,MIUbirc1,[zx;zy]);    
                    COVbirc1 = cat(3,COVbirc1,Qbirc1);    
                end
            end
            for i = 1:NUMpar
                PARrob{r,i}.intenc1 = cat(2,PARrob{r,i}.intenc1,INTENbirc1); PARrob{r,i}.miuc1 = cat(2,PARrob{r,i}.miuc1,MIUbirc1); PARrob{r,i}.covc1 = cat(3,PARrob{r,i}.covc1,COVbirc1);
            end
            
            %% Resample particles
            if MARKresample(1,r) == 1
                PARrobot = PARrob;
                WEIud = PARlikeli / sum(PARlikeli,2);
                PARid = resample(WEIud,NUMpar);
                for i = 1:NUMpar
                    PARrob{r,i}.traj = PARrobot{r,PARid(1,i)}.traj;       PARrob{r,i}.wei = PARrobot{r,PARid(1,i)}.wei;
                    PARrob{r,i}.intenc1 = PARrobot{r,PARid(1,i)}.intenc1; PARrob{r,i}.miuc1 = PARrobot{r,PARid(1,i)}.miuc1; PARrob{r,i}.covc1 = PARrobot{r,PARid(1,i)}.covc1;
                end
            end
        end
%         PARrobrec{1,t} = PARrob;
        if mod(t,200) == 0      % Output for check! 
            t1 = t;
            figure(5);  close;  figure(5);  axis([250 1250 250 1250]);  hold on;    title('Standard PHD-SLAM');
            plot(Landmark_Groundtruth(:,2)',Landmark_Groundtruth(:,3)','ro','MarkerFaceColor','red','MarkerSize',3);
            for r = 1:NUMrobot
                plot(TRAJrec{r,1}(1,1:t1),TRAJrec{r,1}(2,1:t1),'k-');
                plot(Robot_Groundtruth(1,1:t1,1),Robot_Groundtruth(2,1:t1,1),'r-');
                plot(MAPesti{t1,r,tt}(1,:),MAPesti{t1,r,tt}(2,:),'b*');
            end
            pause(0.1);
        end
    end
    TIMErec(1,tt) = toc;
end

% Plot estimated trajectories
figure(1);
hold on;
for r = 1:NUMrobot
    plot(MAPesti{t-1,r,tt}(1,:),MAPesti{t-1,r,tt}(2,:),'b*');
end

% Plot traj estimation error
figure(2);
hold on;
for r = 1:NUMrobot
    h(r) = plot(1:T_obs,PREDerrortraj(r,1:T_obs));
    MEANtraje = mean(TRAJerror,3);
    hh(r) = plot(1:T_obs,MEANtraje(:,r)','--');
end
if NUMrobot > 1
    legend([h(1),hh(1),h(2),hh(2)],'Dead-reckoning-1','Estimated-1','Dead-reckoning-2','Estimated-2');
end

% Plot ospa of map
figure(3);
hold on;
for r = 1:NUMrobot
    MEANospa = mean(OSPAmap,3);
    hh(r) = plot(1:T_obs,MEANospa(:,r)','--');
end




















