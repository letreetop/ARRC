clear all;  close all;  clc;    warning('off');
    
%% System parameters
Pd = 0.9;              % Probability of detection in the field of view
Montecarlo = 1;     
T = 1;                % Sampling time interval
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

%% Record variables
MAPesti = cell(T_obs,NUMrobot,Montecarlo); TRAJrec = cell(NUMrobot,Montecarlo); 
TRAJerror = zeros(T_obs,NUMrobot,Montecarlo); 
TIMErec = zeros(T_obs,Montecarlo);

for tt = 1:Montecarlo
    for r = 1:NUMrobot
        TRAJrec{r,tt} = zeros(3,T_obs); TRAJrec{r,tt}(:,1) = Robot_Groundtruth(:,1,r);
        MAPesti{1,r,tt} = zeros(2,0); 
    end
    for t = 2:T_obs
        for r = 1:NUMrobot
            MAPesti{t,r,tt} = zeros(2,0);
        end
    end
end
pause(0.1);

for tt = 1:Montecarlo
    PARrobrec = cell(1,T_obs);
    %% Generating measurements
    load('Robot_Measurement.mat');
    
    %% Initial SLAM
    MIU = [Robot_Groundtruth(:,1,1)];
    COV = zeros(3,3);  LMcount = zeros(1,0);
    
    %% SLAM begins
    % Attach new landmarks
    for r = 1:NUMrobot
        Zt = Robot_Measurement{r,1};
        xx = MIU( 1+(r-1)*3,1 );     yy = MIU( 2+(r-1)*3,1 );     theta = MIU( 3+(r-1)*3,1 );
        for i = 1:size(Zt,2)
            DIS = Zt(1,i);    ANG = Zt(2,i);
            zx = DIS * cos(theta + ANG) + xx;   zy = DIS * sin(theta + ANG) + yy;       
            figure(1); hold on; plot(zx,zy,'b.','Markersize',8);
            sn = sin(theta + ANG);
            cs = cos(theta + ANG);
            Dgx = [1 0 -Zt(1,i)*sn; 0 1 Zt(1,i)*cs];
            Dgz = [cs, -Zt(1,i)*sn; sn, Zt(1,i)*cs];
            P = COV;
            P = [P(1:3,1:3),     P(1:3,4:end),     P(1:3,1:3)*Dgx';
                 P(4:end,1:3),   P(4:end,4:end),   P(4:end,1:3)*Dgx';
                 Dgx*P(1:3,1:3), Dgx*P(1:3,4:end), Dgx*P(1:3,1:3)*Dgx' + Dgz*R*Dgz'];
            COV = P;    MIU = [MIU',[zx,zy]]';    LMcount = cat(2,LMcount,1);
        end
    end
    
    %% SLAM begins
    MIUall = cell(1,T_obs); COVall = cell(1,T_obs); TIMErec = zeros(1,T_obs);
    for t = 2:T_obs
        tic;
        fprintf('Montecarlo: %d, Time: %d\n',tt,t);
        %% Prediction
        ROBOTstate = MIU( 1:3,1 );  ROBOTcov = COV( 1:3,1:3 );  
        theta = ROBOTstate(3,1);      % Orientation of last time
        ODOrec = Robot_Odometry(:,t,1);
        vel = ODOrec(1,1);    
        vang = ODOrec(2,1);
        F = [T*cos(theta),0;
             T*sin(theta),0;
             0,           T];
        JF = [1, 0, -T*sin(theta)*vel;
              0, 1, T*cos(theta)*vel;
              0, 0, 1];
        JV = [T*cos(theta), 0;
              T*sin(theta), 0;
              0,            T];
        Qtraj = Qodo;
        MIU( 1:3,1 ) = ROBOTstate + F*ODOrec;
        COV( 1:3,1:3 ) = JF * ROBOTcov * JF' + JV * Qtraj * JV';
        
        %% Update
        % Find associated landmarks
        NUMland = (size(MIU,1)-3)/2;    MIUland = zeros(2,NUMland);
        for i = 1:NUMland
            MIUland(:,i) = MIU(3+(i-1)*2+1:3+(i-1)*2+2,1);
        end
        % Solve associations with Murty's algorithm
        ORIrec = cell(1,NUMrobot);  NUMasso = 0;
        for r = 1:NUMrobot
            Zz = Robot_Measurement{r,t};
            ASSOcost = zeros(size(MIUland,2),size(Zz,2));
            for ii = 1:size(MIUland,2)
                for jj = 1:size(Zz,2)
                    xx = MIU((r-1)*3+1,1);  yy = MIU((r-1)*3+2,1);  theta = MIU((r-1)*3+3,1);
                    STATmea = [ Zz(1,jj)*cos(theta+Zz(2,jj)) + xx;Zz(1,jj)*sin(theta+Zz(2,jj)) + yy ];
                    dis = sqrt( (MIUland(:,ii)-STATmea)'*(MIUland(:,ii)-STATmea) );
                    ASSOcost(ii,jj) = dis / 2.5;
                end
            end
            ASSOcost = log(ASSOcost);
            [xxx,yyy] = find(ASSOcost == -Inf);
            if isempty(xxx) ~= 1
                for ii = 1:size(xxx,1)
                    ASSOcost(xxx(ii,1),yyy(ii,1)) = -1000;
                end
            end
            [uasses,nlcost]= mbestwrap_updt_custom(ASSOcost,1);
            ORIrec{1,r} = uasses;   IDnz = find(uasses~=0); NUMasso = NUMasso + size(IDnz,2);
        end
        
        H = zeros( NUMasso*2 , size(MIU,1)  );  RR = zeros(0,0);
        IDnew = cell(4,1);  Zall = zeros(0,1);  Zpred = zeros(0,1); IDcount = 0;
        % Load respectively detected and new landmarks
        for r = 1:NUMrobot
            Zt = Robot_Measurement{r,t};    IDmea = Zt(1,:);    Mt = size(IDmea,2);  IDnew{r,1} = zeros(1,0);  DIFFz = zeros(2,0);
            for i = 1:Mt
                if ismember(i,ORIrec{1,r}) == 1     % If such landmark has been detected
                    Zload = Zt(:,i);  Zload(2,1) = pi_pi(Zload(2,1));  RR = blkdiag(RR,R);
                    Zall = cat(1,Zall,Zload);  POS = find(ORIrec{1,r} == i);   IDcount = IDcount + 1;   LMcount(1,POS) = LMcount(1,POS) + 1;
                    zx = MIU(3+(POS-1)*2+1,1); zy = MIU(3+(POS-1)*2+2,1);
                    xx = MIU( 1+(r-1)*3,1 );     yy = MIU( 2+(r-1)*3,1 );     theta = MIU( 3+(r-1)*3,1 );
                    ANGLEpred = atan2( zy-yy,zx-xx ) - theta;   ANGLEpred = pi_pi(ANGLEpred);  Rpred = sqrt( (zy-yy)^2+(zx-xx)^2 ); 
                    Zpred = cat(1,Zpred,[Rpred;ANGLEpred]);  DIFFz(:,IDcount) = Zload - [Rpred;ANGLEpred];
                    rsqu = (zy-yy)^2+(zx-xx)^2;
                    % Fill in Matrix related to vehicle state
                    H( (IDcount-1)*2+1,(r-1)*3+1 ) = -(zx-xx)/sqrt(rsqu);   H( (IDcount-1)*2+1,(r-1)*3+2 ) = -(zy-yy)/sqrt(rsqu);   H( (IDcount-1)*2+1,(r-1)*3+3 ) = 0;
                    H( (IDcount-1)*2+2,(r-1)*3+1 ) = (zy-yy)/rsqu;          H( (IDcount-1)*2+2,(r-1)*3+2 ) = -(zx-xx)/rsqu;         H( (IDcount-1)*2+2,(r-1)*3+3 ) = -1;
                    % Fill in Matrix related to landmark state
                    H( (IDcount-1)*2+1,3+(POS-1)*2+1 ) = (zx-xx)/sqrt(rsqu);    H( (IDcount-1)*2+1,3+(POS-1)*2+2 ) = (zy-yy)/sqrt(rsqu);
                    H( (IDcount-1)*2+2,3+(POS-1)*2+1 ) = -(zy-yy)/rsqu;         H( (IDcount-1)*2+2,3+(POS-1)*2+2 ) = (zx-xx)/rsqu;
                else    % New landmark
                    IDnew{r,1} = cat(2,IDnew{r,1},i);
                end
            end
        end
        % If there exists detected landmarks, then update
        if IDcount ~= 0    
            S = H * COV * H' + RR;   K = COV * H' * pinv(S);
            COV = (eye(size(COV,2))-K*H) * COV; Zdiff = Zall - Zpred;
            % Invariant estimate
            STATgain = K*Zdiff;  
            GAINref = [STATgain(3,1);STATgain(1:2,1);STATgain(4:end,1)];
            PREref = [MIU(3,1);MIU(1:2,1);MIU(4:end,1)];
            STATEinv = phi(GAINref,PREref);
            MIU = [STATEinv(2:3,1);STATEinv(1,1);STATEinv(4:end)];
        end
        % Attach new landmarks
        for r = 1:NUMrobot
            Zt = Robot_Measurement{r,t}; IDmea = Zt(1,:); NEWid = IDnew{r,1};
            for i = 1:size(NEWid,2)
                xx = MIU( 1+(r-1)*3,1 );     yy = MIU( 2+(r-1)*3,1 );     theta = MIU( 3+(r-1)*3,1 );
                DIS = Zt(1,NEWid(1,i));    ANG = Zt(2,NEWid(1,i));
                zx = DIS * cos(theta + ANG) + xx;   zy = DIS * sin(theta + ANG) + yy;                
                sn = sin(theta + ANG);
                cs = cos(theta + ANG);
                Dgx = [1 0 -Zt(1,NEWid(1,i))*sn; 0 1 Zt(2,NEWid(1,i))*cs];
                Dgz = [cs, -Zt(1,NEWid(1,i))*sn; sn, Zt(1,NEWid(1,i))*cs];
                P = COV;
                P = [P(1:3,1:3),     P(1:3,4:end),     P(1:3,1:3)*Dgx';
                     P(4:end,1:3),   P(4:end,4:end),   P(4:end,1:3)*Dgx';
                     Dgx*P(1:3,1:3), Dgx*P(1:3,4:end), Dgx*P(1:3,1:3)*Dgx' + Dgz*R*Dgz'];
                COV = P;    MIU = [MIU',[zx,zy]]';    LMcount = cat(2,LMcount,1);
            end
        end
        
        % State extraction
        for r = 1:NUMrobot
            % Estimated PHD
            TRAJrec{t,r,tt} = zeros(3,1);
            TRAJrec{t,r,tt} = MIU( (r-1)*3+1:(r-1)*3+3,1 );   % Robot pose record (global coordinate)
            TRAJerror(t,r,tt) = sqrt( (TRAJrec{t,r,tt}(1:2,1) - Robot_Groundtruth(1:2,t,r))' * (TRAJrec{t,r,tt}(1:2,1) - Robot_Groundtruth(1:2,t,r)) );     % Pose estimation error
        end
        
        % Remove undetected landmarks
        if mod(t,4) == 0
            IDremove = find(LMcount <= 2);  IDremove = sort(IDremove,'descend');
            for i = 1:size(IDremove,2)
                POS = IDremove(1,i);
                MIUnew = cat(1,MIU(1:3+(POS-1)*2,1),MIU(3+POS*2+1:end,1));
                COVnew = [COV( 1:3+(POS-1)*2,1:3+(POS-1)*2 ), COV( 1:3+(POS-1)*2,3+POS*2+1:end );
                          COV( 3+POS*2+1:end,1:3+(POS-1)*2 ), COV( 3+POS*2+1:end,3+POS*2+1:end )];
                LMcount(:,POS) = [];
                MIU = MIUnew;   COV = COVnew;
            end
        end
        
        MIUall{1,t} = MIU;  COVall{1,t} = COV;  
        TIMErec(1,t) = toc;
    end
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





















