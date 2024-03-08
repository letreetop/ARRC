function [estimates new] = updateiEKF(estimates,observations,compatibility,H)
% Based on the code of Erik Zamora
% In order to get iEKF-SLAM for the Victoria park dataset,
% Please copy this function to the replace the updateEKF.m file which can be downloaded at https://github.com/ezamorag/EKFSLAM_for_Victoria_Park_Benchmark

global noise
new = [];
STATEpre = estimates.x;
HH = zeros(0,size(compatibility.dH,2)); Zdiff = zeros(0,1); RR = zeros(0,0);
for i=1:size(H,2) 
    j = H(i);
    if j > 0    
        % Update the estimates       
        dH = compatibility.dH(2*j-1:2*j,:);  HH = cat(1,HH,dH);
        innov = [observations.z(1,i)-compatibility.z(1,j); pi_pi(observations.z(2,i)-compatibility.z(2,j))];
        RR = blkdiag(RR,noise.Rz);
        Zdiff = cat(1,Zdiff,innov);
        
        estimates.count(j) = estimates.count(j) + 1;
    else
        % only new features with no neighbours
        if compatibility.AL(i) == 0 
            new = [new i];
        end
    end
end

if size(HH,1) ~= 0
    K =  estimates.P*HH'/(HH* estimates.P*HH' + RR);
    estimates.x = estimates.x + K*Zdiff;
    estimates.P = (eye(3+2*estimates.n)-K*HH)* estimates.P;
end

if size(STATEpre,2) >= 4
    STATgain = estimates.x - STATEpre;  
    GAINref = [STATgain(3,1);STATgain(1:2,1);STATgain(4:end,1)];
    PREref = [STATEpre(3,1);STATEpre(1:2,1);STATEpre(4:end,1)];
    STATEinv = phi(GAINref,PREref);
    estimates.x = [STATEinv(2:3,1);STATEinv(1,1);STATEinv(4:end)];
end



