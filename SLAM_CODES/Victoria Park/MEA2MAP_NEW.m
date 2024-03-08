function [BIRintenall,BIRmiuall,BIRcovall] = MEA2MAP_NEW(NUMrobot,Zpas,POSEesti,PHDmap,t,tt,INTENbir,THRESHOLDmerge,RR,Ptraj)
    BIRintenall = cell(1,NUMrobot);   BIRmiuall = cell(1,NUMrobot);   BIRcovall = cell(1,NUMrobot);
    for r = 1:NUMrobot
        % Initialize landmarks
        BIRinten = zeros(1,0);  BIRmiu = zeros(2,0);  BIRcov = zeros(2,2,0);
        Zpast = Zpas{1,r};
        for i = 1:size(Zpast,2)
            xx = POSEesti(1,r); yy = POSEesti(2,r); theta = POSEesti(3,r);
            LANDmiu = [xx;yy] + [Zpast(1,i)*cos(Zpast(2,i) + theta - pi/2); 
                                 Zpast(1,i)*sin(Zpast(2,i) + theta - pi/2)];
            MARKnew = 1;
            if isempty(PHDmap{t,r,tt}.inten) ~= 1
                for j = 1:size(PHDmap{t,r,tt}.inten,2)      % Check whether there exists similar landmarks
                    dis = sqrt((LANDmiu - PHDmap{t,r,tt}.miu(:,j))'*(LANDmiu - PHDmap{t,r,tt}.miu(:,j)));   
                    if dis <= THRESHOLDmerge
                        MARKnew = 0;
                    end
                end
            end
            sn = sin(Zpast(2,i) + theta - pi/2);
            cs = cos(Zpast(2,i) + theta - pi/2);
            Dgx = [1 0 -Zpast(1,i)*sn; 0 1 Zpast(1,i)*cs];
            Dgz = [cs, -Zpast(1,i)*sn; sn, Zpast(1,i)*cs];
            
            if MARKnew == 1
                BIRinten = cat(2,BIRinten,INTENbir);
                BIRmiu = cat(2,BIRmiu,LANDmiu);
                BIRcov = cat(3,BIRcov,Dgx*Ptraj(:,:,r)*Dgx' + Dgz*RR*Dgz');
            end
        end
        BIRintenall{1,r} = BIRinten;    BIRmiuall{1,r} = BIRmiu;   BIRcovall{1,r} = BIRcov;
    end
end