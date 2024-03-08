function GMPdc1 = PDcompute(ZONEran,FoVANGcen,ZONEang,PDzone,DISpred,ANGpred)

    if DISpred <= ZONEran(1,1) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,1)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,1)/2)  % Zone1
        GMPdc1 = PDzone(1,1);
    elseif DISpred <= ZONEran(1,2) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,2)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,2)/2)  % Zone2
        GMPdc1 = PDzone(1,2);
    elseif DISpred <= ZONEran(1,3) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,3)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,3)/2)  % Zone3
        GMPdc1 = PDzone(1,3);
    elseif DISpred <= ZONEran(1,4) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,4)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,4)/2)  % Zone4
        GMPdc1 = PDzone(1,4);
    elseif DISpred <= ZONEran(1,5) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,5)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,5)/2)  % Zone5
        GMPdc1 = PDzone(1,5);
    elseif DISpred <= ZONEran(1,6) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,6)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,6)/2)  % Zone6
        GMPdc1 = PDzone(1,6);
    elseif DISpred <= ZONEran(1,7) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,7)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,7)/2)  % Zone7
        GMPdc1 = PDzone(1,7);
    elseif DISpred <= ZONEran(1,8) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,8)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,8)/2)  % Zone8
        GMPdc1 = PDzone(1,8);
    elseif DISpred <= ZONEran(1,9) && ANGpred*180/pi >= (FoVANGcen - ZONEang(1,9)/2) && ANGpred*180/pi <= (FoVANGcen + ZONEang(1,9)/2)  % Zone9
        GMPdc1 = PDzone(1,9);
    else
        GMPdc1 = PDzone(1,10);
    end

end