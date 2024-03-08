function Figureoutput(t,tt,IDfig,GPSdata,TRAJrec,MAPesti,NUMrobot)

    load('TRAJref');   % load('LANDMARKref');
    figure(IDfig);
    % subplot(2,2,3);
    axis([-200 200 -150 350]);
    % axis([-100 0 -60 0])
    hold on;
%     GPSn = GPSdata(:,1:t);  IDo = find(GPSn(1,:) == 0); GPSn(:,IDo) = [];
    GPSn = GPSdata;  IDo = find(GPSn(1,:) == 0); GPSn(:,IDo) = [];
    plot(GPSn(1,:),GPSn(2,:),'ro','MarkerSize',1,'MarkerFaceColor','r');
%     plot(TRAJref(1,1:t),TRAJref(2,1:t),'co','MarkerSize',2,'MarkerFaceColor','c');

    hold on;
    % plot(TRAJdr(1,1:T_obs),TRAJdr(2,1:T_obs),'bo','MarkerSize',2,'MarkerFaceColor','blue');

    for r = 1:NUMrobot
        TRAJr = TRAJrec{r,tt};
        plot(TRAJr(1,1:t),TRAJr(2,1:t),'k.','MarkerSize',4,'MarkerFaceColor','black');
        plot(MAPesti{t,r,tt}(1,:),MAPesti{t,r,tt}(2,:),'b.','MarkerSize',8);
    end
%     figure(IDfig);
%     hold on;
%     plot(LANDMARKref(1,:),LANDMARKref(2,:),'bo');
end