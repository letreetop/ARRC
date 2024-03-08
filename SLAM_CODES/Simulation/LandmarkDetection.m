function [EMPTrec,LANDMARKdetpt,LANDMARKdet,MAPnumtrue] = LandmarkDetection(NUMrobot,T_obs,Robot_Groundtruth,Landmark_Groundtruth,RANGminfov,RANGmaxfov,ANGLEminfov,ANGLEmaxfov)

    LANDMARKdet = cell(NUMrobot,T_obs);
    LANDMARKdetpt = cell(NUMrobot,T_obs);
    EMPTrec = zeros(1,0);
    for t = 1:T_obs
        for r = 1:NUMrobot
            xx = Robot_Groundtruth(1,t,r);  yy = Robot_Groundtruth(2,t,r);  theta = Robot_Groundtruth(3,t,r);   theta = pi_pi(theta);
            LMtem = zeros(1,0);
            for i = 1:size(Landmark_Groundtruth,1)
                zx = Landmark_Groundtruth(i,2); zy = Landmark_Groundtruth(i,3);
                DIS = sqrt( (xx-zx)^2 + (yy-zy)^2 );
                ANG = atan2( zy-yy,zx-xx ) - theta; ANG = pi_pi(ANG);
                ANG = ANG * 180 / pi;
                if DIS >= RANGminfov(1,r) && DIS <= RANGmaxfov(1,r) && ANG >= ANGLEminfov(1,r) && ANG <= ANGLEmaxfov(1,r)
                    LMtem = cat(2,LMtem,Landmark_Groundtruth(i,1));
                end
            end
            if isempty(LMtem) == 1
                EMPTrec = cat(2,EMPTrec,t);
            end
            LANDMARKdetpt{r,t} = LMtem;
            LANDMARKdet{r,t} = LANDMARKdetpt{r,t};
            if t ~= 1
                LANDMARKdet{r,t} = cat(2,LANDMARKdet{r,t},LANDMARKdet{r,t-1});
                LANDMARKdet{r,t} = unique(LANDMARKdet{r,t});
                LANDMARKdet{r,t} = sort(LANDMARKdet{r,t});
            end
        end
    end

    MAPnumtrue = zeros(1,T_obs);
    for ts = 1:T_obs
        DETMAP = [];
        for r = 1:NUMrobot
            DETMAP = [DETMAP,LANDMARKdet{r,ts}];
        end
        DETMAP = unique(DETMAP);
        MAPnumtrue(1,ts) = size(DETMAP,2);
    end

end