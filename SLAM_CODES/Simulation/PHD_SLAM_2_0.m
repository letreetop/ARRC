clear all;  close all;  clc;    warning('off');
    
%% System parameters
Pdo = 0.95;
Pd = Pdo;              % Probability of detection in the field of view
Ps = 0.95;
Montecarlo = 1;     
T = 1;                % Sampling time interval
NUMpar = 5;            % Number of particles
THRESHOLDprun = 1e-5;   % Prun threshold
THRESHOLDmerge = 2.5;   % Merge threshold
THRESHOLDdetect = 0.5;  % Landmark detecting threshold
INTENbir = 0.05;         % Intensity of birth GMs
NUMclutter = 5;
MARKunifibirth = 1;     % If use unifi birth, set to 1; otherwise set to 0

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
            PARrob{r,i}.P = diag([0.001,0.001,0.0001*pi/180]).^2;   PARrob{r,i}.wei = 1/NUMpar;        % weight of particle
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
        TRAJallpar = zeros(3,NUMpar);
        for r = 1:NUMrobot
            Qtraj = Qodo(:,:,r);
            ZM = Robot_Measurement{r,t};    
            PARlikeli = zeros(1,NUMpar); % Likelihood just computed based on class 1
            for i = 1:NUMpar
                PARrob{r,i}.P = diag([0.001,0.001,0.0001*pi/180]).^2;
                %% Generating birth landmarks
                if MARKunifibirth ~= 1
                    xx = PARrob{r,i}.traj(1,1);     yy = PARrob{r,i}.traj(2,1);     theta = PARrob{r,i}.traj(3,1);      
                    ZP = Robot_Measurement{r,t-1};
                    INTENbirc1 = zeros(1,0);    MIUbirc1 = zeros(2,0);  COVbirc1 = zeros(2,2,0);
                    for j = 1:size(ZP,2)
                        DIS = ZP(1,j);    ANG = ZP(2,j);
                        zx = DIS * cos(theta + ANG) + xx;   zy = DIS * sin(theta + ANG) + yy;
                        rsqu = (zy-yy)^2+(zx-xx)^2;
                        HH = [ (zx-xx)/sqrt(rsqu),(zy-yy)/sqrt(rsqu);
                               -(zy-yy)/rsqu,     (zx-xx)/rsqu];
                        Qbir = pinv( HH' * pinv(R) * HH );
                        MIUtot = PARrob{r,i}.miuc1;
                        MATRIXdis = repmat([zx;zy],1,size(MIUtot,2)) - MIUtot;
                        DISall = sqrt( MATRIXdis(1,:).^2+MATRIXdis(2,:).^2 );
                        if min(DISall) >= 5     % If the measurement is not close to any existing landmark/target
                            INTENbirc1 = cat(2,INTENbirc1,INTENbir);
                            MIUbirc1 = cat(2,MIUbirc1,[zx;zy]);    
                            COVbirc1 = cat(3,COVbirc1,Qbir);    
                        end
                    end
                    PARrob{r,i}.intenc1 = cat(2,PARrob{r,i}.intenc1,INTENbirc1); PARrob{r,i}.miuc1 = cat(2,PARrob{r,i}.miuc1,MIUbirc1); PARrob{r,i}.covc1 = cat(3,PARrob{r,i}.covc1,COVbirc1);
                end
                %% Prediction for trajectory
                theta = PARrob{r,i}.traj(3,1);      % Orientation of last time
                ODOrec = Robot_Odometry(:,t,r);   vel = ODOrec(1,1);    vang = ODOrec(2,1);     % Odometry at this time
                F = [T*cos(theta),0;
                     T*sin(theta),0;
                     0,           T];
                JF = [1, 0, -T*sin(theta)*vel;
                      0, 1, T*cos(theta)*vel;
                      0, 0, 1];
                JV = [T*cos(theta), 0;
                      T*sin(theta), 0;
                      0,            T];
                PREDodo = PARrob{r,i}.traj + F*ODOrec;      PARrob{r,i}.traj = PREDodo;     % Trajectory prediction strictly follows the motion model
                PARrob{r,i}.P = JF * PARrob{r,i}.P * JF' + JV * Qtraj * JV';                % Propagation of covariance
                PREDtraj = PARrob{r,i}.traj;    PREDP = PARrob{r,i}.P;
                
                xx = PARrob{r,i}.traj(1,1);     yy = PARrob{r,i}.traj(2,1);     theta = PARrob{r,i}.traj(3,1);      PP = PARrob{r,i}.P;
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
                NUMgm = size(INTENc1,2);    INTENaug = INTENc1;
                MIUaug = cat(1,MIUc1,repmat([xx;yy;theta],1,NUMgm));    COVaug = zeros(5,5,NUMgm);
                for j = 1:size(COVc1,3)      
                    COVaug(:,:,j) = blkdiag( COVc1(:,:,j),PP );
                end
                
                
                %% Update
                if isempty(IDc1fov) ~= 1 || size(ZM,2) == 0           % Update carry out if and only if there are static landmarks inside the current FoV
                    %% STEP 1: Importance sampling
                    % Step 1A: Static map selection
                    Zpre = zeros(2,NUMgm);  GMPD = zeros(1,NUMgm);  GMid = zeros(1,0);
                    HH = zeros(2,5,NUMgm);      SS = zeros(2,2,NUMgm);  KK = zeros(5,2,NUMgm);  COVuptem = zeros(5,5,NUMgm);
                    for j = 1:NUMgm
                        zx = MIUc1(1,j); zy = MIUc1(2,j); 
                        ANGpred = atan2( zy-yy,zx-xx ) - theta;  ANGpred = pi_pi(ANGpred);  DISpred = sqrt( (zy-yy)^2+(zx-xx)^2 );     Zpre(:,j) = [DISpred;ANGpred];
                        rsqu = (zy-yy)^2+(zx-xx)^2;
                        HH(:,:,j) = [ (zx-xx)/sqrt(rsqu),(zy-yy)/sqrt(rsqu),-(zx-xx)/sqrt(rsqu),-(zy-yy)/sqrt(rsqu),0;
                                      -(zy-yy)/rsqu,     (zx-xx)/rsqu,      (zy-yy)/rsqu,       -(zx-xx)/rsqu,      -1];
                        SS(:,:,j) = HH(:,:,j) * COVaug(:,:,j) * HH(:,:,j)' + R;  KK(:,:,j) = COVaug(:,:,j) * HH(:,:,j)' * pinv(SS(:,:,j));
                        COVuptem(:,:,j) = (eye(5)-KK(:,:,j)*HH(:,:,j)) * COVaug(:,:,j);     % Used only for record
                        if DISpred <= RANGmaxfov(1,r) && ANGpred*180/pi >= (ANGLEminfov(1,r)) && ANGpred*180/pi <= (ANGLEmaxfov(1,r))   % In the current FoV
                            GMid = cat(2,GMid,j);
                            PDtem = PDcompute(ZONEran(r,:),FoVANGcen(1,r),ZONEang(r,:),PDzone,DISpred,ANGpred);    GMPD(1,j) = PDtem;
                        end
                    end
                    % Step 1B: Update static map
                    LIKELIall = zeros(1,size(ZM,2));        % Used to normalize PHD
                    INTENup = zeros(1,0);               MIUup = zeros(5,0);     COVup = zeros(5,5,0);   
                    for j = 1:size(ZM,2)
                        LIKELIgm = zeros(1,NUMgm);  zreal = ZM(1:2,j);
                        for k = 1:NUMgm
                            zpred = Zpre(:,k); zdiff = zreal-zpred; zdiff(2,1) = pi_pi(zdiff(2,1));
                            LIKELIgm(1,k) = GMPD(1,k) * (det(2*pi*SS(:,:,k)))^(-1) * exp( -0.5*zdiff'*pinv(SS(:,:,k))*zdiff ) * INTENaug(1,k) + 1e-299;
                            INTENup = cat(2,INTENup,LIKELIgm(1,k));
                            COVup = cat(3,COVup,COVuptem(:,:,k));
                            MIUup = cat(2,MIUup, MIUaug(:,k) + KK(:,:,k)*zdiff );
                        end
                        LIKELIall(1,j) = sum(LIKELIgm,2);
                    end
                    LIKELIcomb = sum(LIKELIall,1);   
                    % Step 1C: Compute the final PHD
                    NUMgm = size(INTENc1,2);
                    for j = 1:size(ZM,2)
                        Id = (j-1)*NUMgm+1:j*NUMgm;
                        INTENup(1,Id) = INTENup(1,Id) ./(INTENclutter+LIKELIcomb(1,j));
                    end
                    
                    % Step 1D: Separate landmark and trajectory state
                    WEItraj = INTENup / sum(INTENup,2); TRAJup = MIUup(3:5,:);      COVtraj = zeros(3,3,size(COVup,3)); 
                    for j = 1:size(COVup,3)
                        COVtraj(:,:,j) = COVup(3:5,3:5,j); 
                    end
                    IDw = find(INTENup >= 0.5);         
                    
                    % Discarding first-part PHD
                    WEIt = WEItraj(1,IDw) / sum(WEItraj(1,IDw),2);    
                    TRAJup = TRAJup(:,IDw); COVtraj = COVtraj(:,:,IDw); 
                    % ----------
                    
                    OMEGAtraj = zeros(3,3,size(IDw,2)); qup = zeros(3,size(IDw,2));
                    for j = 1:size(IDw,2)
                        OMEGAtraj(:,:,j) = WEIt(1,j) * pinv(COVtraj(:,:,j));
                        qup(:,j) = WEIt(1,j) * pinv(COVtraj(:,:,j)) * TRAJup(:,j);
                    end
                    % Estimate trajectory and covariance matrix
                    TRAJesti = sum(repmat(WEIt,3,1).*TRAJup,2);  Pup = zeros(3,3);
                    for j = 1:size(qup,2)
                        Pup = Pup + WEIt(1,j) * (COVtraj(:,:,j) + (TRAJesti-TRAJup(:,j))*(TRAJesti-TRAJup(:,j))');
                    end

                    if isempty(IDw) ~= 1
                        MARKss = 1;
                        while (MARKss) 
                            PARrob{r,i}.P = Pup;    PARrob{r,i}.traj = TRAJesti + chol(Pup)*randn(3,1);    TRAJsam = PARrob{r,i}.traj;  
                            WEIimportance = 1/sqrt(det(2*pi*PREDP)) * exp( -0.5*(TRAJsam-PREDtraj)'*pinv(PREDP)*(TRAJsam-PREDtraj) ) ...
                                          / ( 1/sqrt(det(2*pi*Pup)) * exp( -0.5*(TRAJsam-TRAJesti)'*pinv(Pup)*(TRAJsam-TRAJesti) ) );
                            if isnan(WEIimportance) == 1
                                MARKss = 1;
                            else
                                MARKss = 0;
                            end
                        end
                                  
                        xx = PARrob{r,i}.traj(1,1); yy = PARrob{r,i}.traj(2,1); theta = PARrob{r,i}.traj(3,1);      % Revised vehicle trajectory information
                    else
                        WEIimportance = 1;
                    end
                    TRAJallpar(:,i) = PARrob{r,i}.traj;
                    Pd = Pdo;
                    %% STEP 2: Update based on the newly sampled trajectory
                     % Step 2 - Pre: Find the Gaussian inside the FoV
                    MIU = PARrob{r,i}.miuc1;  IDc1fov = zeros(1,0);  
                    for j = 1:size(MIU,2)     % Find landmarks that inside the FoV and then attach with Pd
                        zx = MIU(1,j); zy = MIU(2,j);
                        ANGpred = atan2( zy-yy,zx-xx ) - theta;  ANGpred = pi_pi(ANGpred);  DISpred = sqrt( (zy-yy)^2+(zx-xx)^2 ); 
                        if DISpred <= RANGmaxfov(1,r)+0.6 && ANGpred*180/pi >= (ANGLEminfov(1,r))-3 && ANGpred*180/pi <= (ANGLEmaxfov(1,r))+3   % In the current FoV
                            IDc1fov = cat(2,IDc1fov,j);    
                        end
                    end
                    IDc1out = setdiff([1:size(MIU,2)],IDc1fov);     % Static landmarks that outside the FoV, which not be involved in the current iteration
                    INTENc1out = PARrob{r,i}.intenc1(1,IDc1out);    MIUc1out = PARrob{r,i}.miuc1(:,IDc1out);    COVc1out = PARrob{r,i}.covc1(:,:,IDc1out);
                    % Landmarks that remain static
                    INTENc1 = PARrob{r,i}.intenc1(1,IDc1fov);     MIUc1 = PARrob{r,i}.miuc1(:,IDc1fov);   COVc1 = PARrob{r,i}.covc1(:,:,IDc1fov);
                    NUMgm = size(INTENc1,2);    INTENaug = INTENc1;
                    MIUaug = cat(1,MIUc1,repmat([xx;yy;theta],1,NUMgm));    COVaug = zeros(5,5,NUMgm);
                    for j = 1:size(COVc1,3)       % Do we really need to propagate the covariance matrices of landmarks?
                        COVaug(:,:,j) = blkdiag( COVc1(:,:,j),PP );
                    end
                    % Step 2A: Static map selection
                    Zpre = zeros(2,NUMgm);  GMPD = zeros(1,NUMgm);  GMid = zeros(1,0);
                    HH = zeros(2,2,NUMgm);      SS = zeros(2,2,NUMgm);  KK = zeros(2,2,NUMgm);  COVuptem = zeros(2,2,NUMgm);
                    for j = 1:NUMgm
                        zx = MIUc1(1,j); zy = MIUc1(2,j); 
                        ANGpred = atan2( zy-yy,zx-xx ) - theta;  ANGpred = pi_pi(ANGpred);  DISpred = sqrt( (zy-yy)^2+(zx-xx)^2 );     Zpre(:,j) = [DISpred;ANGpred];
                        rsqu = (zy-yy)^2+(zx-xx)^2;
                        HH(:,:,j) = [ (zx-xx)/sqrt(rsqu),(zy-yy)/sqrt(rsqu);
                                      -(zy-yy)/rsqu,     (zx-xx)/rsqu];
                        SS(:,:,j) = HH(:,:,j) * COVc1(:,:,j) * HH(:,:,j)' + R;  KK(:,:,j) = COVc1(:,:,j) * HH(:,:,j)' * pinv(SS(:,:,j));
                        COVuptem(:,:,j) = (eye(2)-KK(:,:,j)*HH(:,:,j)) * COVc1(:,:,j);     % Used only for record
                        if DISpred <= RANGmaxfov(1,r) && ANGpred*180/pi >= (ANGLEminfov(1,r)) && ANGpred*180/pi <= (ANGLEmaxfov(1,r))   % In the current FoV
                            GMid = cat(2,GMid,j);
                            PDtem = PDcompute(ZONEran(r,:),FoVANGcen(1,r),ZONEang(r,:),PDzone,DISpred,ANGpred);    GMPD(1,j) = PDtem;
                        end
                    end
                    % Step 2B: Compute particle likelihood
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
                    PARlikeli(1,i) = prod(LIKELIpz,2) * WEIimportance + 1e-99;
                    % Step 2C: Update static map
                    LIKELIall = zeros(1,size(ZM,2));        % Used to normalize PHD
                    INTENup = zeros(1,0);               MIUup = zeros(2,0);     COVup = zeros(2,2,0);   
                    INTENudc1 = (1-GMPD).*INTENc1;      MIUudc1 = MIUc1;        COVudc1 = COVc1;
                    for j = 1:size(ZM,2)
                        LIKELIgm = zeros(1,NUMgm);  zreal = ZM(1:2,j);
                        for k = 1:NUMgm
                            zpred = Zpre(:,k); zdiff = zreal-zpred; zdiff(2,1) = pi_pi(zdiff(2,1));
                            LIKELIgm(1,k) = GMPD(1,k) * (det(2*pi*SS(:,:,k)))^(-1) * exp( -0.5*zdiff'*pinv(SS(:,:,k))*zdiff ) * INTENc1(1,k) + 1e-299;
                            INTENup = cat(2,INTENup,LIKELIgm(1,k));
                            COVup = cat(3,COVup,COVuptem(:,:,k));
                            MIUup = cat(2,MIUup, MIUc1(:,k) + KK(:,:,k)*zdiff );
                        end
                        LIKELIall(1,j) = sum(LIKELIgm,2);
                    end
                    LIKELIcomb = sum(LIKELIall,1);   
                    % Step 2D: Compute the final PHD
                    NUMgm = size(INTENc1,2);
                    for j = 1:size(ZM,2)
                        Id = (j-1)*NUMgm+1:j*NUMgm;
                        INTENup(1,Id) = INTENup(1,Id) ./(INTENclutter+LIKELIcomb(1,j));
                    end
                    
                    INTENupc1 = min(INTENup,0.995);     MIUupc1 = MIUup(1:2,:);     COVupc1 = zeros(2,2,size(COVup,3));
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
            
            %% Birth PHD generate
            if MARKunifibirth == 1
                ZM = Robot_Measurement{r,t};
                INTENbirc1 = zeros(1,0);    MIUbirc1 = zeros(2,0);  COVbirc1 = zeros(2,2,0);
                for j = 1:size(ZM,2)
                    xx = POSEesti(1,r);     yy = POSEesti(2,r);     theta = POSEesti(3,r);
                    DIS = ZM(1,j);    ANG = ZM(2,j);
                    zx = DIS * cos(theta + ANG) + xx;   zy = DIS * sin(theta + ANG) + yy;
                    rsqu = (zy-yy)^2+(zx-xx)^2;
                    HH = [ (zx-xx)/sqrt(rsqu),(zy-yy)/sqrt(rsqu);
                           -(zy-yy)/rsqu,     (zx-xx)/rsqu];
                    Qbir = pinv( HH' * pinv(R) * HH );
                    MIUtot = PARrob{r,MAXid(1,1)}.miuc1;
                    MATRIXdis = repmat([zx;zy],1,size(MIUtot,2)) - MIUtot;
                    DISall = sqrt( MATRIXdis(1,:).^2+MATRIXdis(2,:).^2 );
                    if min(DISall) >= 5     % If the measurement is not close to any existing landmark/target
                        INTENbirc1 = cat(2,INTENbirc1,INTENbir);
                        MIUbirc1 = cat(2,MIUbirc1,[zx;zy]);    
                        COVbirc1 = cat(3,COVbirc1,Qbir);    
                    end
                end
                for i = 1:NUMpar
                    PARrob{r,i}.intenc1 = cat(2,PARrob{r,i}.intenc1,INTENbirc1); PARrob{r,i}.miuc1 = cat(2,PARrob{r,i}.miuc1,MIUbirc1); PARrob{r,i}.covc1 = cat(3,PARrob{r,i}.covc1,COVbirc1);
                end
            end
            
            %% Resample particles
            if MARKresample(1,r) == 1
                PARrobot = PARrob;
                WEIud = PARlikeli / sum(PARlikeli,2);   IDp = find(WEIud == max(WEIud));
                PARid = resample(WEIud,NUMpar);     PARid = cat(2,IDp(1,1),PARid);  PARid = PARid(1,1:NUMpar);
                for i = 1:NUMpar
                    PARrob{r,i}.traj = PARrobot{r,PARid(1,i)}.traj;  PARrob{r,i}.traj(3,1) = pi_pi(PARrob{r,i}.traj(3,1));
                    PARrob{r,i}.P = PARrobot{r,PARid(1,i)}.P;
                    PARrob{r,i}.wei = 1/NUMpar;
                    PARrob{r,i}.intenc1 = PARrobot{r,PARid(1,i)}.intenc1; PARrob{r,i}.miuc1 = PARrobot{r,PARid(1,i)}.miuc1; PARrob{r,i}.covc1 = PARrobot{r,PARid(1,i)}.covc1;
                end
            end
        end
%         PARrobrec{1,t} = PARrob;
        if mod(t,200) == 0      % Output for check! 
            t1 = t;
            figure(5);  close;  figure(5);  axis([250 1250 250 1250]);  hold on;    title('Proposed method');
            plot(Landmark_Groundtruth(:,2)',Landmark_Groundtruth(:,3)','ro','MarkerFaceColor','red','MarkerSize',3);
            for r = 1:NUMrobot
                plot(TRAJrec{r,tt}(1,1:t1),TRAJrec{r,tt}(2,1:t1),'k-');
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




















