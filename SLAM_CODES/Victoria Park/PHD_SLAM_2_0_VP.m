clear all;close all;clc;warning('off');

%% System parameters
Pd = 0.9;              % Probability of detection in the field of view
MARKrec = 0;             % Decide whether to record data (it will take a lot of memeory space!)
Montecarlo = 1;     
T_obs = 61945;          % Observation time
NUMpar = 1;             % Number of particles
THRESHOLDprun = 1e-5;   % 8e-10;   % Prun threshold
THRESHOLDmerge = 2;   % Merge threshold
THRESHOLDdetect = 0.5;  % Landmark detecting threshold
MAXGMnum = 100;
INTENbir = 0.01;         % Intensity of birth GMs
NUMclutter = 0.2;
THRESHOLDpd = 5e-1;    % Threshold for extracting the Gaussian components to build the proposal density 
R = diag([1 3*pi/180]).^2;
PARrem = min(5,NUMpar);
NUMrobot = 1;
INTENclutter = 1.8394e-06;
UNIFIbirth = 1;         % 1: All particles use same birth; 0: Each particle generate its own birth
DISper = 70;

%% Vehicle parameters
VEHL = 2.83;    VEHa = 3.78;    VEHb = 0.5;     VEH = 0.76;
Qmodel = diag([2,2,0.015]);
POSini = [-67.6493; -41.7142; 36*pi/180]; %Initial condition

% General FoV
ANGLEmaxfov = 180*ones(1,5);
ANGLEminfov = 0*ones(1,5);
RANGminfov = 0*ones(1,5);
RANGmaxfov = 75*ones(1,5);
FoVANGcen = 0.5*(ANGLEminfov+ANGLEmaxfov);
ZONEang = zeros(5,10);  ZONEran = zeros(5,10);
for r = 1:NUMrobot
    for i = 1:10
        ZONEang(r,i) = ANGLEminfov(1,r) + (ANGLEmaxfov(1,r)-ANGLEminfov(1,r))/10*i;
        ZONEran(r,i) = RANGminfov(1,r) + (RANGmaxfov(1,r)-RANGminfov(1,r))/10*i;
    end
end
PDzone = [Pd,Pd,Pd*0.88,Pd*0.88,Pd*0.88,Pd*0.32,Pd*0.32,Pd*0.32,Pd*0.1,Pd*0.01];

%% Record variables
MAPesti = cell(T_obs,NUMrobot,Montecarlo);
TRAJrec = cell(NUMrobot,Montecarlo);
PHDmap = cell(T_obs,NUMrobot,Montecarlo);
OSPAmap = zeros(T_obs,NUMrobot,Montecarlo);
TRAJerror = zeros(T_obs,NUMrobot,Montecarlo);
NUMlandmark = zeros(T_obs,NUMrobot,Montecarlo);
for tt = 1:Montecarlo
    for r = 1:NUMrobot
        TRAJrec{r,tt} = zeros(3,T_obs);
        TRAJrec{r,tt}(:,1) = POSini;
        MAPesti{1,r,tt} = zeros(2,0);
    end
    for t = 1:T_obs
        for r = 1:NUMrobot
            MAPesti{t,r,tt} = zeros(2,0);
            PHDmap{t,r,tt}.inten = zeros(1,0);  PHDmap{t,r,tt}.miu = zeros(2,0);    PHDmap{t,r,tt}.cov = zeros(2,2,0);
        end
    end
end

load('VPdatareshape.mat');
figure(1);
axis([-200 200 -150 350]);
hold on;
GPSn = GPSdata(:,1:T_obs);  IDo = find(GPSn(1,:) == 0); GPSn(:,IDo) = [];
plot(GPSn(1,:),GPSn(2,:),'ro','MarkerSize',2,'MarkerFaceColor','r');

%% SLAM start
for tt = 1:Montecarlo
    load('MEAori.mat');     % Load Dataset
    %% Filter initial
    PARrob = cell(NUMpar,NUMrobot);
    POSEesti = zeros(3,NUMrobot);  
    for i = 1:NUMpar
        for r = 1:NUMrobot
            PARrob{i,r}.traj = POSini;    PARrob{i,r}.wei = 1/NUMpar;    PARrob{i,r}.P = diag([0.01,0.01,0.001*pi/180]).^2;  
            PARrob{i,r}.inten = zeros(1,0);   PARrob{i,r}.miu = zeros(2,0);  PARrob{i,r}.cov = zeros(2,2,0); 
        end
    end
    Zpast = cell(1,NUMrobot);    
    for r = 1:NUMrobot
        POSEesti(:,r) = POSini;
        Zpast{1,r} = Z{r,1};
    end
    BIRintenrec = cell(T_obs,NUMrobot);     BIRmiurec = cell(T_obs,NUMrobot);     BIRcovrec = cell(T_obs,NUMrobot);
    Prr = repmat(diag([0.01,0.01,0.001*pi/180]).^2,1,1,NUMrobot);
    [BIRintenall,BIRmiuall,BIRcovall] = MEA2MAP_NEW(NUMrobot,Zpast,POSEesti,PHDmap,1,1,INTENbir,THRESHOLDmerge,R,Prr);
%     figure(1); hold on; plot(BIRmiuall{1,1}(1,:),BIRmiuall{1,1}(2,:),'bo');
    for r = 1:NUMrobot
        BIRintenrec{1,r} = BIRintenall{1,r};   BIRmiurec{1,r} = BIRmiuall{1,r};   BIRcovrec{1,r} = BIRcovall{1,r};
        for i = 1:NUMpar
            PARrob{i,r}.inten = BIRintenall{1,r};  PARrob{i,r}.miu = BIRmiuall{1,r};  PARrob{i,r}.cov = BIRcovall{1,r};
        end
    end
    
    %% SLAM begin
    IDfig = 2;  
    PARrobrec = cell(1,T_obs);
    
    TIMEexe = zeros(1,T_obs);
    for t = 2:T_obs
        tic;
        fprintf('Montecarlo: %d, Time: %d\n',tt,t);
        % Prediction
        PARlikeliall = zeros(NUMrobot,NUMpar);
        for r = 1:NUMrobot
            ZM = Z{r,t}; 
            % Find the measurements
            PARlikeli = zeros(1,NUMpar);
            for i = 1:NUMpar
                PARrob{i,r}.P = diag([0.01,0.01,0.001*pi/180]).^2;  
                % Prediction of trajectory
                ucur = u(:,t-1);     vt = ucur(1,1);    at = ucur(2,1);     Vc = vt / (1-tan(at)*VEH/VEHL);
                th = PARrob{i,r}.traj(3,1);
                tmp1 = VEHa*sin(th) + VEHb*cos(th);
                tmp2 = VEHa*cos(th) - VEHb*sin(th);
                
                Dfx = [1 0 -T*Vc*(sin(th) + 1/VEHL*tan(at)*(VEHa*cos(th) - VEHb*sin(th))); 
                       0 1  T*Vc*(cos(th) - 1/VEHL*tan(at)*(VEHa*sin(th) + VEHb*cos(th)));
                       0 0  1];
                
                PARrob{i,r}.traj = PARrob{i,r}.traj + [T*Vc*(cos(th) - 1/VEHL*tan(at)*tmp1); 
                                                       T*Vc*(sin(th) + 1/VEHL*tan(at)*tmp2);
                                                       T*Vc/VEHL*tan(at)];
                PARrob{i,r}.P = Dfx * PARrob{i,r}.P * Dfx' + Qmodel;
                
                % Prediction of landmarks (static, just copy!)
                if t ~= 2
                    % Attach new landmarks
                    PARrob{i,r}.inten = cat(2,PARrob{i,r}.inten,BIRintenrec{t-1,r});
                    PARrob{i,r}.miu = cat(2,PARrob{i,r}.miu,BIRmiurec{t-1,r});
                    PARrob{i,r}.cov = cat(3,PARrob{i,r}.cov,BIRcovrec{t-1,r});
                else
                    PARrob{i,r}.inten = PARrob{i,r}.inten;
                    PARrob{i,r}.miu = PARrob{i,r}.miu;
                    PARrob{i,r}.cov = PARrob{i,r}.cov;
                end
                xx = PARrob{i,r}.traj(1,1);   yy = PARrob{i,r}.traj(2,1);   theta = minimizedAngle(PARrob{i,r}.traj(3,1));
                % Update
                % Find the landmarks in current FoV
                NUMgm = size(PARrob{i,r}.inten,2);  
                if NUMgm ~= 0 && isempty(ZM) ~= 1
                    %% Important sampling
                    GMPd = zeros(1,NUMgm);
                    GMid = zeros(1,0);  Zpre = zeros(2,NUMgm);  HH = zeros(2,5,0);  SS = zeros(2,2,0);  KK = zeros(5,2,0);
                    COVupt = zeros(5,5,0);                    
                    % Compute the likelihood of the particle                    
                    LIKELIpz = zeros(1,size(ZM,2));     INTENpredsum = 1;
                    INTENup = zeros(1,0);               MIUup = zeros(5,0);      COVup = zeros(5,5,0);
                    for j = 1:size(ZM,2)
                        LIKELIpf = zeros(1,0);  LIKELIgm = zeros(1,0);
                        for k = 1:NUMgm
                            zx = PARrob{i,r}.miu(1,k);  zy = PARrob{i,r}.miu(2,k);
                            if j == 1   % Check whether it is inside the FoV
                                ANGLEpred = minimizedAngle(atan2( zy-yy,zx-xx ) - theta + pi/2);
                                rrpred = sqrt( (zy-yy)^2+(zx-xx)^2 );
                                zdiff = ZM(1:2,j) - [rrpred;ANGLEpred];  zdiff(2,1) = minimizedAngle(zdiff(2,1));   Zpre(:,k) = [rrpred;ANGLEpred];
                                if rrpred <= RANGmaxfov(1,r)+5 && ANGLEpred*180/pi >= (ANGLEminfov(1,r)-10) && ANGLEpred*180/pi <= (ANGLEmaxfov(1,r)+10)   % In the current FoV
                                    GMid = cat(2,GMid,k);                                    
                                    rsqu = (zy-yy)^2+(zx-xx)^2;
                                    H = [ (zx-xx)/sqrt(rsqu), (zy-yy)/sqrt(rsqu), -(zx-xx)/sqrt(rsqu), -(zy-yy)/sqrt(rsqu), 0;
                                           -(zy-yy)/rsqu,     (zx-xx)/rsqu,       (zy-yy)/rsqu,        -(zx-xx)/rsqu,       -1 ];
                                    COVaug = blkdiag(PARrob{i,r}.cov(:,:,k),PARrob{i,r}.P);
                                    S = H * COVaug * H' + R;  K = COVaug * H' * pinv(S);
                                    COVt = (eye(5)-K*H)*COVaug;
                                    HH = cat(3,HH,H);	SS = cat(3,SS,S);   KK = cat(3,KK,K);
                                    COVupt = cat(3,COVupt,COVt);
                                    if rrpred <= ZONEran(r,1) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,1)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,1)/2)  % Zone1
                                        GMPd(1,k) = PDzone(1,1);
                                    elseif rrpred <= ZONEran(r,2) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,2)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,2)/2)  % Zone2
                                        GMPd(1,k) = PDzone(1,2);
                                    elseif rrpred <= ZONEran(r,3) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,3)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,3)/2)  % Zone3
                                        GMPd(1,k) = PDzone(1,3);
                                    elseif rrpred <= ZONEran(r,4) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,4)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,4)/2)  % Zone4
                                        GMPd(1,k) = PDzone(1,4);
                                    elseif rrpred <= ZONEran(r,5) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,5)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,5)/2)  % Zone5
                                        GMPd(1,k) = PDzone(1,5);
                                    elseif rrpred <= ZONEran(r,6) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,6)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,6)/2)  % Zone6
                                        GMPd(1,k) = PDzone(1,6);
                                    elseif rrpred <= ZONEran(r,7) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,7)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,7)/2)  % Zone7
                                        GMPd(1,k) = PDzone(1,7);
                                    elseif rrpred <= ZONEran(r,8) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,8)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,8)/2)  % Zone8
                                        GMPd(1,k) = PDzone(1,8);
                                    elseif rrpred <= ZONEran(r,9) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,9)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,9)/2)  % Zone9
                                        GMPd(1,k) = PDzone(1,9);
                                    else
                                        GMPd(1,k) = PDzone(1,10);
                                    end
                                end
                            end
                            if ismember(k,GMid) == 1    % If it is inside the FoV
                                zdiff = ZM(1:2,j) - Zpre(:,k);  zdiff(2,1) = minimizedAngle(zdiff(2,1));
                                POSid = find( GMid == k );
                                if abs(zdiff(1,1)) <= 1 && abs(zdiff(2,1)) <= (3*pi/180)
                                    COVup = cat(3,COVup,COVupt(:,:,POSid));
                                    MIUup = cat(2,MIUup, [zx;zy;xx;yy;theta]+KK(:,:,POSid)*zdiff );
                                    LIKELIgmtem = GMPd(1,GMid(1,POSid)) * (det(2*pi*SS(:,:,POSid)))^(-0.5) * exp( -0.5*zdiff'*pinv(SS(:,:,POSid))*zdiff ) * PARrob{i,r}.inten(1,k) + 1e-299;
                                    LIKELIgm = cat(2,LIKELIgm,LIKELIgmtem);
                                end
                            end
                        end
                        LIKELIall = sum(LIKELIgm,2);
                        INTENup = cat( 2,INTENup,min(LIKELIgm./(INTENclutter+LIKELIall),0.999) );
                    end
                    
                    % Convert the augmented state back to separated sapces
                    MIUupaug = MIUup;   COVupaug = COVup;   
                    MIUtraj = MIUupaug(3:5,:);  COVtraj = zeros(3,3,size(COVupaug,3));
                    for j = 1:size(COVupaug,3)
                        COVtraj(:,:,j) = COVupaug(3:5,3:5,j);
                    end
                    if isempty(MIUup) ~= 1
                        IDw = find(INTENup >= THRESHOLDpd); MIUtraj = MIUtraj(:,IDw);   COVtraj = COVtraj(:,:,IDw);
                        WEIt = INTENup(1,IDw) / sum(INTENup(1,IDw),2);
                        TRAJestiCU = sum(repmat(WEIt,3,1).*MIUtraj,2);  PupCU = zeros(3,3); qup = zeros(3,size(IDw,2));
                        OMEGAup = zeros(3,3,size(IDw,2));
                        for j = 1:size(WEIt,2)
                            PupCU = PupCU + WEIt(1,j) * (COVtraj(:,:,j) + (TRAJestiCU-MIUtraj(:,j))*(TRAJestiCU-MIUtraj(:,j))');
                            OMEGAup(:,:,j) = WEIt(1,j) * pinv(COVtraj(:,:,j));
                            qup(:,j) = WEIt(1,j) * pinv(COVtraj(:,:,j)) * MIUtraj(:,j);
                        end
                        PupCI = pinv(sum(OMEGAup(:,:,j),3));  TRAJestiCI = PupCI * sum(qup,2);
                        if isempty(IDw) ~= 1
                            if i == 1
                                PARrob{i,r}.traj = TRAJestiCU;        PARrob{i,r}.P = PupCU;
                            else
                                PARrob{i,r}.traj = TRAJestiCU + chol(PupCU)/100*randn(3,1);        
                                PARrob{i,r}.P = PupCU;
                            end
                        end
                    end
                    if isempty(IDw) ~= 1
                        xx = PARrob{i,r}.traj(1,1); yy = PARrob{i,r}.traj(2,1); theta = PARrob{i,r}.traj(3,1);      % Revise trajectory information
                    end
                    
                    %% Real update!
                    GMPd = zeros(1,NUMgm);
                    GMid = zeros(1,0);  Zpre = zeros(2,NUMgm);  HH = zeros(2,2,0);  SS = zeros(2,2,0);  KK = zeros(2,2,0);
                    COVupt = zeros(2,2,0);                    
                    % Compute the likelihood of the particle                    
                    LIKELIpz = zeros(1,size(ZM,2));     INTENpredsum = 1;
                    INTENup = zeros(1,0);               MIUup = zeros(2,0);      COVup = zeros(2,2,0);
                    for j = 1:size(ZM,2)
                        LIKELIpf = zeros(1,0);  LIKELIgm = zeros(1,0);
                        for k = 1:NUMgm
                            zx = PARrob{i,r}.miu(1,k);  zy = PARrob{i,r}.miu(2,k);
                            if j == 1   % Check whether it is inside the FoV
                                ANGLEpred = minimizedAngle(atan2( zy-yy,zx-xx ) - theta + pi/2);
                                rrpred = sqrt( (zy-yy)^2+(zx-xx)^2 );
                                zdiff = ZM(1:2,j) - [rrpred;ANGLEpred];  zdiff(2,1) = minimizedAngle(zdiff(2,1));   Zpre(:,k) = [rrpred;ANGLEpred];
                                if rrpred <= RANGmaxfov(1,r)+5 && ANGLEpred*180/pi >= (ANGLEminfov(1,r)-10) && ANGLEpred*180/pi <= (ANGLEmaxfov(1,r)+10)   % In the current FoV
                                    GMid = cat(2,GMid,k);                                    
                                    rsqu = (zy-yy)^2+(zx-xx)^2;
                                    H = [ (zx-xx)/sqrt(rsqu), (zy-yy)/sqrt(rsqu);
                                           -(zy-yy)/rsqu,     (zx-xx)/rsqu];
                                    COVaug = PARrob{i,r}.cov(:,:,k);
                                    S = H * COVaug * H' + R;  K = COVaug * H' * pinv(S);
                                    COVt = (eye(2)-K*H)*COVaug;
                                    HH = cat(3,HH,H);	SS = cat(3,SS,S);   KK = cat(3,KK,K);
                                    COVupt = cat(3,COVupt,COVt);
                                    if rrpred <= ZONEran(r,1) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,1)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,1)/2)  % Zone1
                                        GMPd(1,k) = PDzone(1,1);
                                    elseif rrpred <= ZONEran(r,2) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,2)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,2)/2)  % Zone2
                                        GMPd(1,k) = PDzone(1,2);
                                    elseif rrpred <= ZONEran(r,3) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,3)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,3)/2)  % Zone3
                                        GMPd(1,k) = PDzone(1,3);
                                    elseif rrpred <= ZONEran(r,4) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,4)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,4)/2)  % Zone4
                                        GMPd(1,k) = PDzone(1,4);
                                    elseif rrpred <= ZONEran(r,5) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,5)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,5)/2)  % Zone5
                                        GMPd(1,k) = PDzone(1,5);
                                    elseif rrpred <= ZONEran(r,6) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,6)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,6)/2)  % Zone6
                                        GMPd(1,k) = PDzone(1,6);
                                    elseif rrpred <= ZONEran(r,7) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,7)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,7)/2)  % Zone7
                                        GMPd(1,k) = PDzone(1,7);
                                    elseif rrpred <= ZONEran(r,8) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,8)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,8)/2)  % Zone8
                                        GMPd(1,k) = PDzone(1,8);
                                    elseif rrpred <= ZONEran(r,9) && ANGLEpred*180/pi >= (FoVANGcen(1,r) - ZONEang(r,9)/2) && ANGLEpred*180/pi <= (FoVANGcen(1,r) + ZONEang(r,9)/2)  % Zone9
                                        GMPd(1,k) = PDzone(1,9);
                                    else
                                        GMPd(1,k) = PDzone(1,10);
                                    end
                                end
                            end
                            if ismember(k,GMid) == 1    % If it is inside the FoV
                                zdiff = ZM(1:2,j) - Zpre(:,k);  zdiff(2,1) = minimizedAngle(zdiff(2,1));
                                POSid = find( GMid == k );
                                if rrpred <= DISper && ZM(1,j) <= DISper
                                    LIKELIpf = cat( 2,LIKELIpf,GMPd(1,k)/sqrt(det(2*pi*R))*exp( -0.5*(zdiff)'*pinv(R)*zdiff ) * PARrob{i,r}.inten(1,k) );
                                end
                                if abs(zdiff(1,1)) <= 2 && abs(zdiff(2,1)) <= (3*pi/180)
                                    COVup = cat(3,COVup,COVupt(:,:,POSid));
                                    MIUup = cat(2,MIUup, [zx;zy]+KK(:,:,POSid)*zdiff );
                                    LIKELIgmtem = GMPd(1,GMid(1,POSid)) * (det(2*pi*SS(:,:,POSid)))^(-0.5) * exp( -0.5*zdiff'*pinv(SS(:,:,POSid))*zdiff ) * PARrob{i,r}.inten(1,k) + 1e-299;
                                    LIKELIgm = cat(2,LIKELIgm,LIKELIgmtem);
                                end
                            end
                        end
                        LIKELIpz(1,j) = INTENclutter + sum(LIKELIpf,2);
                        LIKELIall = sum(LIKELIgm,2);
                        INTENup = cat( 2,INTENup,min(LIKELIgm./(INTENclutter+LIKELIall),0.999) );
                    end
                    PARlikeli(1,i) = prod(LIKELIpz,2)*INTENpredsum + 1e-99;
                    
                    GMidout = setdiff(1:NUMgm,GMid);
                    % Remain PHDs outside current FoV
                    INTENout = PARrob{i,r}.inten(1,GMidout);  MIUout = PARrob{i,r}.miu(:,GMidout);  COVout = PARrob{i,r}.cov(:,:,GMidout);
                    % Compute undetected PHD
                    INTENud = PARrob{i,r}.inten(1,GMid) .* (1-GMPd(1,GMid));   MIUud = PARrob{i,r}.miu(:,GMid);           COVud = PARrob{i,r}.cov(:,:,GMid);
                    INTENup = cat(2,INTENud,INTENup);   MIUup = cat(2,MIUud,MIUup);     COVup = cat(3,COVud,COVup);
                    [w1,x1,P1] = gaus_prune(INTENup,MIUup,COVup,THRESHOLDprun);  [w2,x2,P2] = gaus_merge(w1',x1,P1,THRESHOLDmerge);  
                    WEImerge = w2'; MIUmerge = x2;  COVmerge = P2;
                    WEImerge = min(WEImerge,0.999);
                    INTENall = cat(2,INTENout,WEImerge);  MIUall = cat(2,MIUout,MIUmerge);  COVall = cat(3,COVout,COVmerge);
                    PARrob{i,r}.inten = INTENall;     PARrob{i,r}.miu = MIUall;    PARrob{i,r}.cov = COVall;
                else
                    PARlikeli(1,i) = 1/NUMpar;
                end
            end
            PARlikeliall(r,:) = PARlikeli;
        end
         %% Find local estimated PHD and resample particle IDs
        PHDmaplocal = cell(1,NUMrobot); IDres = cell(1,NUMrobot);
        for r = 1:NUMrobot
            % Estimated PHD
            PARlikeli = PARlikeliall(r,:);
            MAXid = find(PARlikeli == max(PARlikeli));
            PHDmaplocal{1,r}.inten = PARrob{MAXid(1,1),r}.inten; PHDmaplocal{1,r}.miu = PARrob{MAXid(1,1),r}.miu; PHDmaplocal{1,r}.cov = PARrob{MAXid(1,1),r}.cov;
            % Resample IDs
            WEIud = PARlikeli / sum(PARlikeli,2);
            % Standard resample
            IDres{1,r} = resample(WEIud,NUMpar);
            IDres{1,r} = cat( 2,MAXid(1,1),IDres{1,r}(1,1:NUMpar-1) );      % Ensure the particle with largest weight is included
        end
         %% State estimation
        PPesti = zeros(3,3,NUMrobot);
        for r = 1:NUMrobot
            % Trajectory estimation
            PARlikeli = PARlikeliall(r,:);
            MAXid = find(PARlikeli == max(PARlikeli));
            POSEesti(:,r) = PARrob{MAXid(1,1),r}.traj;         % Estimated pose (local coordinate)
            TRAJrec{r,tt}(:,t) = zeros(3,1);
            TRAJrec{r,tt}(1:2,t) = PARrob{MAXid(1,1),r}.traj(1:2,1);   % Robot pose record (global coordinate)
            TRAJrec{r,tt}(3,t) = PARrob{MAXid(1,1),r}.traj(3,1);
            PPesti = PARrob{MAXid(1,1),r}.P;
            % Landmark estimation
            IDmap = find(PARrob{MAXid(1,1),r}.inten >= THRESHOLDdetect);  
            PHDmap{t,r,tt}.inten = PARrob{MAXid(1,1),r}.inten;
            PHDmap{t,r,tt}.miu = PARrob{MAXid(1,1),r}.miu;
            PHDmap{t,r,tt}.cov = PARrob{MAXid(1,1),r}.cov;
            NUMlandmark(t,r,tt) = size(IDmap,2);    % Number of detected landmarks
            MAPesti{t,r,tt} = PARrob{MAXid(1,1),r}.miu(:,IDmap);    
            if MARKrec == 0
                PHDmap{t-1,r,tt}.inten = zeros(1,0); PHDmap{t-1,r,tt}.miu = zeros(2,0); PHDmap{t-1,r,tt}.cov = zeros(2,2,0);
                MAPesti{t-1,r,tt} = zeros(2,0);
            end
        end
        % Generate birth for next time
        if UNIFIbirth == 1
            Zpast = cell(1,NUMrobot);    
            for r = 1:NUMrobot
                Zpast{1,r} = Z{r,t};
            end
            [BIRintenall,BIRmiuall,BIRcovall] = MEA2MAP_NEW(NUMrobot,Zpast,POSEesti,PHDmap,t,tt,INTENbir,THRESHOLDmerge,R,PPesti);
            for r = 1:NUMrobot
                BIRintenrec{t,r} = BIRintenall{1,r};   BIRmiurec{t,r} = BIRmiuall{1,r};   BIRcovrec{t,r} = BIRcovall{1,r};
            end
        end
            
        %% Resample particles
        if isempty(Z{r,t}) ~= 1
            PARrobot = PARrob;   
            for r = 1:NUMrobot
                PARid = IDres{1,r};
                for i = 1:NUMpar
                    PARrob{i,r}.traj = PARrobot{PARid(1,i),r}.traj; 
                    PARrob{i,r}.wei = 1/NUMpar;  PARrob{i,r}.Q = PARrobot{PARid(1,i),r}.P;
                    PARrob{i,r}.inten = PARrobot{PARid(1,i),r}.inten; PARrob{i,r}.miu = PARrobot{PARid(1,i),r}.miu; PARrob{i,r}.cov = PARrobot{PARid(1,i),r}.cov;
                end
            end
        end
        
        TIMEexe(1,t) = toc;
        % Output
        if t == 100 || mod(t,1000) == 0
            IDfig = tt + 1;
            figure(IDfig);  clf;
            Figureoutput(t,tt,IDfig,GPSdata,TRAJrec,MAPesti,NUMrobot);
            pause(0.1);
        end
        if MARKrec == 1
            PARrobrec{1,t} = PARrob;
        end
    end
end


figure(1);
hold on;
for r = 1:NUMrobot
    TRAJr = TRAJrec{r,1};
    plot(TRAJr(1,:),TRAJr(2,:),'ko','MarkerSize',2,'MarkerFaceColor','black');
    if r == 1
        plot(MAPesti{T_obs,r,tt}(1,:),MAPesti{T_obs,r,tt}(2,:),'b*');
    else
        plot(MAPesti{T_obs,r,tt}(1,:),MAPesti{T_obs,r,tt}(2,:),'g*');
    end
end














