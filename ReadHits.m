%  Hits file (gateHits.dat)
%
% Each line is a hit and the columns represent :
%
%     Column 1 : ID of the run (i.e. time-slice)
%     Column 2 : ID of the event
%     Column 3 : ID of the primary particle whose descendant generated this hit
%     Column 4 : ID of the source which emitted the primary particle
%     Columns 5 to N+4: Volume IDs at each level of the hierarchy of a system,
%       so the number of columns depends on the system used.
%
% For cylindricalPET system N=6 :
%
%     Column 5 : ID of volume attached to the "base" level of the system
%     Column 6 : ID of volume attached to the "rsector" level of the system
%     Column 7 : ID of volume attached to the "module" level of the system
%     Column 8 : ID of volume attached to the "submodule" level of the system
%     Column 9 : ID of volume attached to the "crystal" level of the system
%     Column 10 : ID of volume attached to the "layer" level of the system
%
% For SPECTHead system N=3 :
%
%     Column 5 : ID of volume attached to the "base" level of the system
%     Column 6 : ID of volume attached to the "crystal" level of the system
%     Column 7 : ID of volume attached to the "pixel" level of the system
%     Column N+5 : Time stamp of the hit
%     Column N+6 : Energy deposited by the hit (which may be given as a percentage of the initial particle)
%     Column N+7 : Range of particle which has generated the hit
%     Column N+8, N+9 ,N+10 : XYZ position of the hit in the world referential
%     Column N+11 : Geant4 code of the particle which has generated the hit (11 for Electrons & 22 for Photons)
%     Column N+12 : ID of the particle which has generated the hit
%     Column N+13 : ID of the mother of the particle which has generated the hit
%     Column N+14 : ID of the photon giving the particle which has generated the hit
%     Column N+15 : Number of Compton interactions in phantoms before reaching the detector
%     Column N+16 : Number of Rayleigh interactions in phantoms before reaching the detector
%     Column N+17 : Name of the process which has generated the hit
%     Column N+18 : Name of the last volume where a Compton effect occurred
%     Column N+19 : Name of the last volume where a Rayleigh effect occurred
%
% For the next sections, the system will be set as a cylindricalPET system,
%   so that the number of lines concerning the Volume ID of each level will always be 5.
%
EventPointerFirst=zeros(1,1);
EventPointerLast=zeros(1,1);
EventIDByHit=zeros(20000000,1);
LastEvent=0;
%
fid=fopen('AsciiFileHits.dat','r');
nHits=0;
tline='Temp';
while ischar(tline)
    %disp(tline)
    tline = fgetl(fid);
    nHits=nHits+1;
    ThisHit=nHits;
    if (mod(nHits,10000)==0)
        disp(['nHits=',num2str(nHits)]);
    end
    if (numel(tline)>1)
        [A,count,errmsg,nextindex]=sscanf(tline,'%d %d %d %d %d %d %d');
        ThisEvent=A(2)+1;
        EventIDByHit(nHits)=ThisEvent;
        if (ThisEvent>LastEvent)
            EventPointerFirst(ThisEvent)=ThisHit;
            if (LastEvent>0)
                EventPointerLast(LastEvent)=LastHit;
            end
        end
        LastEvent=ThisEvent;
        LastHit=nHits;
    end
    if (nHits>20000000)
        break
    end
end
EventIDByHit=EventIDByHit(1:nHits);
EventPointerLast(LastEvent)=LastHit;
nEvents=A(2)+1;
disp(['nEvents=',num2str(nEvents),' nHits=',num2str(nHits)]);
fclose(fid);
%
% Check Hits
% CheckMean=zeros(nEvents,1);
% CheckStd=zeros(nEvents,1);
% for iEvent=1:nEvents
%     if (EventPointerFirst(iEvent)>0)
%         Check=EventIDByHit(EventPointerFirst(iEvent):EventPointerLast(iEvent));
%         CheckMean(iEvent)=mean(Check);
%         CheckStd(iEvent)=std(Check);
%         disp(['iEvent=',num2str(iEvent),' Check Mean=',num2str(mean(Check)),' Std=',num2str(std(Check))]);
%     end
% end
%
MaxHits=max(EventPointerLast-EventPointerFirst+1);
%
%%%
%
Trans.nTrans=zeros(nEvents,1);
Absorption.nAbsorption=zeros(nEvents,1);
Compton.nCompton=zeros(nEvents,1);
Cerenkov.nCerenkov=zeros(nEvents,1);
Phot.nPhot=zeros(nEvents,1);
Msc.nMsc=zeros(nEvents,1);
Rayl.nRayl=zeros(nEvents,1);
eIoni.neIoni=zeros(nEvents,1);
eBrem.neBrem=zeros(nEvents,1);
Scintillation.nScintillation=zeros(nEvents,1);
%
TotalEnergy=zeros(nEvents,1);
NEventComptons=zeros(nEvents,1);
XFirstHit=zeros(nEvents,2);
XPhotocapture=zeros(nEvents,1);
YPhotocapture=zeros(nEvents,1);
ZPhotocapture=zeros(nEvents,1);
%
fid=fopen('AsciiFileHits.dat','r');
%
nHit=0;
for iEvent=1:nEvents
    nHitsThisEvent=EventPointerLast(iEvent)-EventPointerFirst(iEvent)+1;
    if (EventPointerFirst(iEvent)>0)
        nTrans=0;
        nAbsorption=0;
        nCompton=0;
        nPhot=0;
        nMsc=0;
        nRayl=0;
        neIoni=0;
        neBrem=0;
        nCerenkov=0;
        nScintillation=0;
        %
        for iHit=1:nHitsThisEvent
            nHit=nHit+1;
%             if (iHit==1)
%                 disp(['Event ',num2str(iEvent),' firstHit=',num2str(nHit), ...
%                     ' EventPointerFirst=',num2str(EventPointerFirst(iEvent))]);
%             end
%             if (iHit==nHitsThisEvent)
%                 disp(['Event ',num2str(iEvent),' LastHit=',num2str(nHit),' EventPointerLast=',num2str(EventPointerLast(iEvent))]);
%             end
            %
            tline = fgetl(fid);
            if (mod(nHit,10000)==0)
                disp(['nHit=',num2str(nHit)]);
            end
            if (numel(tline)>1)
                test=strfind(tline,char(32));
                tline(test)=' ';
                [A,count,errmsg,nextindex]=sscanf(tline,'%d %d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d %d %d');
                tleft=tline(nextindex:end);
                WhiteSpaceNum=find(isstrprop(tleft,'wspace'));
                Str1=sscanf(tleft(1:(WhiteSpaceNum(1)-1)),'%s');
                Str2=sscanf(tleft((WhiteSpaceNum(1)+1):WhiteSpaceNum(2)-1),'%s');
                Str3=sscanf(tleft((WhiteSpaceNum(2)+1):end),'%s');
                %
                switch Str1
                    case 'Transportation'
                        nTrans=nTrans+1;
                        Trans.nTrans(iEvent)=nTrans;
                        %
                    case 'OpticalAbsorption'
                        nAbsorption=nAbsorption+1;
                        Absorption.nAbsorption(iEvent)=nAbsorption;
                    case 'compt'
                        nCompton=nCompton+1;
                        Compton.nCompton(iEvent)=nCompton;
                    case 'phot'
                        nPhot=nPhot+1;
                        Phot.nPhot(iEvent)=nPhot;
                    case 'Cerenkov'
                        nCerenkov=nCerenkov+1;
                        Cerenkov.nCerenkov(iEvent)=nCerenkov;
                    case 'msc'
                        nMsc=nMsc+1;
                        Msc.nMsc(iEvent)=nMsc;
                    case 'Rayl'
                        nRayl=nRayl+1;
                        Rayl.nRayl(iEvent)=nRayl;
                    case 'eIoni'
                        neIoni=neIoni+1;
                        eIoni.neIoni(iEvent)=neIoni;
                    case 'eBrem'
                        neBrem=neBrem+1;
                        eBrem.neBrem(iEvent)=neBrem;
                    case 'Scintillation'
                        nScintillation=nScintillation+1;
                        Scintillation.nScintillation(iEvent)=nScintillation;
                    otherwise
                        disp(['Unusual Str1 Content: ',Str1]);
                        disp(tline);
                        %                    disp('here');
                end
            end
        end
    end
end
fclose(fid);
%
Total.nTrans=nTrans;
Total.nAbsorption=nAbsorption;
Total.nCompton=nCompton;
Total.nPhot=nPhot;
Total.nCerenkov=nCerenkov;
Total.nMsc=nMsc;
Total.nRayl=nRayl;
Total.neIoni=neIoni;
Total.neBram=neBrem;
Total.nScintillation=nScintillation;
%
%%%
%
nTransMax=max(Trans.nTrans);
Trans.nHit=zeros(nEvents,nTransMax);
Trans.RunID=zeros(nEvents,nTransMax);
Trans.EventID=zeros(nEvents,nTransMax);
Trans.PrimaryParticleID=zeros(nEvents,nTransMax);
Trans.SourceID=zeros(nEvents,nTransMax);
Trans.BaseID=zeros(nEvents,nTransMax);
Trans.CrystalID=zeros(nEvents,nTransMax);
Trans.PixelID=zeros(nEvents,nTransMax);
Trans.Time=zeros(nEvents,nTransMax);
Trans.Energy=zeros(nEvents,nTransMax);
Trans.Range=zeros(nEvents,nTransMax);
Trans.XPos=zeros(nEvents,nTransMax);
Trans.YPos=zeros(nEvents,nTransMax);
Trans.ZPos=zeros(nEvents,nTransMax);
Trans.GParticle=zeros(nEvents,nTransMax);
Trans.HitParticleID=zeros(nEvents,nTransMax);
Trans.MotherParticleID=zeros(nEvents,nTransMax);
Trans.PhotonID=zeros(nEvents,nTransMax);
Trans.NAbsorption=zeros(nEvents,nTransMax);
Trans.NRayleigh=zeros(nEvents,nTransMax);
%
nAbsorptionMax=max(Absorption.nAbsorption);
Absorption.nHit=zeros(nEvents,nAbsorptionMax);
Absorption.RunID=zeros(nEvents,nAbsorptionMax);
Absorption.EventID=zeros(nEvents,nAbsorptionMax);
Absorption.PrimaryParticleID=zeros(nEvents,nAbsorptionMax);
Absorption.SourceID=zeros(nEvents,nAbsorptionMax);
Absorption.BaseID=zeros(nEvents,nAbsorptionMax);
Absorption.CrystalID=zeros(nEvents,nAbsorptionMax);
Absorption.PixelID=zeros(nEvents,nAbsorptionMax);
Absorption.Time=zeros(nEvents,nAbsorptionMax);
Absorption.Energy=zeros(nEvents,nAbsorptionMax);
Absorption.Range=zeros(nEvents,nAbsorptionMax);
Absorption.XPos=zeros(nEvents,nAbsorptionMax);
Absorption.YPos=zeros(nEvents,nAbsorptionMax);
Absorption.ZPos=zeros(nEvents,nAbsorptionMax);
Absorption.GParticle=zeros(nEvents,nAbsorptionMax);
Absorption.HitParticleID=zeros(nEvents,nAbsorptionMax);
Absorption.MotherParticleID=zeros(nEvents,nAbsorptionMax);
Absorption.PhotonID=zeros(nEvents,nAbsorptionMax);
Absorption.NAbsorption=zeros(nEvents,nAbsorptionMax);
Absorption.NRayleigh=zeros(nEvents,nAbsorptionMax);
%
nComptonMax=max(Compton.nCompton);
Compton.nHit=zeros(nEvents,nComptonMax);
Compton.RunID=zeros(nEvents,nComptonMax);
Compton.EventID=zeros(nEvents,nComptonMax);
Compton.PrimaryParticleID=zeros(nEvents,nComptonMax);
Compton.SourceID=zeros(nEvents,nComptonMax);
Compton.BaseID=zeros(nEvents,nComptonMax);
Compton.CrystalID=zeros(nEvents,nComptonMax);
Compton.PixelID=zeros(nEvents,nComptonMax);
Compton.Time=zeros(nEvents,nComptonMax);
Compton.Energy=zeros(nEvents,nComptonMax);
Compton.Range=zeros(nEvents,nComptonMax);
Compton.XPos=zeros(nEvents,nComptonMax);
Compton.YPos=zeros(nEvents,nComptonMax);
Compton.ZPos=zeros(nEvents,nComptonMax);
Compton.GParticle=zeros(nEvents,nComptonMax);
Compton.HitParticleID=zeros(nEvents,nComptonMax);
Compton.MotherParticleID=zeros(nEvents,nComptonMax);
Compton.PhotonID=zeros(nEvents,nComptonMax);
Compton.NCompton=zeros(nEvents,nComptonMax);
Compton.NRayleigh=zeros(nEvents,nComptonMax);
%
nCerenkovMax=max(Cerenkov.nCerenkov);
Cerenkov.nHit=zeros(nEvents,nCerenkovMax);
Cerenkov.RunID=zeros(nEvents,nCerenkovMax);
Cerenkov.EventID=zeros(nEvents,nCerenkovMax);
Cerenkov.PrimaryParticleID=zeros(nEvents,nCerenkovMax);
Cerenkov.SourceID=zeros(nEvents,nCerenkovMax);
Cerenkov.BaseID=zeros(nEvents,nCerenkovMax);
Cerenkov.CrystalID=zeros(nEvents,nCerenkovMax);
Cerenkov.PixelID=zeros(nEvents,nCerenkovMax);
Cerenkov.Time=zeros(nEvents,nCerenkovMax);
Cerenkov.Energy=zeros(nEvents,nCerenkovMax);
Cerenkov.Range=zeros(nEvents,nCerenkovMax);
Cerenkov.XPos=zeros(nEvents,nCerenkovMax);
Cerenkov.YPos=zeros(nEvents,nCerenkovMax);
Cerenkov.ZPos=zeros(nEvents,nCerenkovMax);
Cerenkov.GParticle=zeros(nEvents,nCerenkovMax);
Cerenkov.HitParticleID=zeros(nEvents,nCerenkovMax);
Cerenkov.MotherParticleID=zeros(nEvents,nCerenkovMax);
Cerenkov.PhotonID=zeros(nEvents,nCerenkovMax);
Cerenkov.NCerenkov=zeros(nEvents,nCerenkovMax);
Cerenkov.NRayleigh=zeros(nEvents,nCerenkovMax);
%
nPhotMax=max(Phot.nPhot);
Phot.nHit=zeros(nEvents,nPhotMax);
Phot.RunID=zeros(nEvents,nPhotMax);
Phot.EventID=zeros(nEvents,nPhotMax);
Phot.PrimaryParticleID=zeros(nEvents,nPhotMax);
Phot.SourceID=zeros(nEvents,nPhotMax);
Phot.BaseID=zeros(nEvents,nPhotMax);
Phot.CrystalID=zeros(nEvents,nPhotMax);
Phot.PixelID=zeros(nEvents,nPhotMax);
Phot.Time=zeros(nEvents,nPhotMax);
Phot.Energy=zeros(nEvents,nPhotMax);
Phot.Range=zeros(nEvents,nPhotMax);
Phot.XPos=zeros(nEvents,nPhotMax);
Phot.YPos=zeros(nEvents,nPhotMax);
Phot.ZPos=zeros(nEvents,nPhotMax);
Phot.GParticle=zeros(nEvents,nPhotMax);
Phot.HitParticleID=zeros(nEvents,nPhotMax);
Phot.MotherParticleID=zeros(nEvents,nPhotMax);
Phot.PhotonID=zeros(nEvents,nPhotMax);
Phot.NPhot=zeros(nEvents,nPhotMax);
Phot.NRayleigh=zeros(nEvents,nPhotMax);
%
fid=fopen('AsciiFileHits.dat','r');
%
nHit=0;
for iEvent=1:nEvents
    nHitsThisEvent=EventPointerLast(iEvent)-EventPointerFirst(iEvent)+1;
    if (EventPointerFirst(iEvent)>0)
        nTrans=0;
        nAbsorption=0;
        nCompton=0;
        nPhot=0;
        nMsc=0;
        nRayl=0;
        neIoni=0;
        neBrem=0;
        nCerenkov=0;
        nScintillation=0;
        %
        for iHit=1:nHitsThisEvent
            nHit=nHit+1;
%             disp(['nHit=',num2str(nHit)]);
            %
            tline = fgetl(fid);
            if (mod(nHit,10000)==0)
                disp(['nHit=',num2str(nHit)]);
            end
            if (numel(tline)>1)
                test=strfind(tline,char(32));
                tline(test)=' ';
                [A,count,errmsg,nextindex]=sscanf(tline,'%d %d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d %d %d');
                tleft=tline(nextindex:end);
                WhiteSpaceNum=find(isstrprop(tleft,'wspace'));
                Str1=sscanf(tleft(1:(WhiteSpaceNum(1)-1)),'%s');
                Str2=sscanf(tleft((WhiteSpaceNum(1)+1):WhiteSpaceNum(2)-1),'%s');
                Str3=sscanf(tleft((WhiteSpaceNum(2)+1):end),'%s');
                %
                %
                switch Str1
                    case 'Transportation'
                        nTrans=nTrans+1;
                        Trans.nHit(iEvent,nTrans)=nHit;
                        Trans.RunID(iEvent,nTrans)=A(1);
%                         if (iEvent==2)
%                             disp('here');
%                         end
                        Trans.EventID(iEvent,nTrans)=A(2);
                        Trans.PrimaryParticleID(iEvent,nTrans)=A(3);
                        Trans.SourceID(iEvent,nTrans)=A(4);
                        Trans.BaseID(iEvent,nTrans)=A(5);
                        Trans.CrystalID(iEvent,nTrans)=A(6);
                        Trans.PixelID(iEvent,nTrans)=A(7);
                        Trans.Time(iEvent,nTrans)=A(8);
                        Trans.Energy(iEvent,nTrans)=A(9);
                        Trans.Range(iEvent,nTrans)=A(10);
                        Trans.XPos(iEvent,nTrans)=A(11);
                        Trans.YPos(iEvent,nTrans)=A(12);
                        Trans.ZPos(iEvent,nTrans)=A(13);
                        Trans.GParticle(iEvent,nTrans)=A(14);
                        Trans.HitParticleID(iEvent,nTrans)=A(15);
                        Trans.MotherParticleID(iEvent,nTrans)=A(16);
                        Trans.PhotonID(iEvent,nTrans)=A(17);
                        Trans.NCompton(iEvent,nTrans)=A(18);
                        Trans.NRayleigh(iEvent,nTrans)=A(19);
                    case 'OpticalAbsorption'
                        nAbsorption=nAbsorption+1;
                        Absorption.nHit(iEvent,nAbsorption)=nHit;
                        Absorption.RunID(iEvent,nAbsorption)=A(1);
                        Absorption.EventID(iEvent,nAbsorption)=A(2);
                        Absorption.PrimaryParticleID(iEvent,nAbsorption)=A(3);
                        Absorption.SourceID(iEvent,nAbsorption)=A(4);
                        Absorption.BaseID(iEvent,nAbsorption)=A(5);
                        Absorption.CrystalID(iEvent,nAbsorption)=A(6);
                        Absorption.PixelID(iEvent,nAbsorption)=A(7);
                        Absorption.Time(iEvent,nAbsorption)=A(8);
                        Absorption.Energy(iEvent,nAbsorption)=A(9);
                        Absorption.Range(iEvent,nAbsorption)=A(10);
                        Absorption.XPos(iEvent,nAbsorption)=A(11);
                        Absorption.YPos(iEvent,nAbsorption)=A(12);
                        Absorption.ZPos(iEvent,nAbsorption)=A(13);
                        Absorption.GParticle(iEvent,nAbsorption)=A(14);
                        Absorption.HitParticleID(iEvent,nAbsorption)=A(15);
                        Absorption.MotherParticleID(iEvent,nAbsorption)=A(16);
                        Absorption.PhotonID(iEvent,nAbsorption)=A(17);
                        Absorption.NAbsorption(iEvent,nAbsorption)=A(18);
                        Absorption.NRayleigh(iEvent,nAbsorption)=A(19);
                    case 'compt'
                        nCompton=nCompton+1;
                        Compton.nCompton(iEvent)=nCompton;
                        Compton.nHit(iEvent,nCompton)=nHit;
                        Compton.RunID(iEvent,nCompton)=A(1);
                        Compton.EventID(iEvent,nCompton)=A(2);
                        Compton.PrimaryParticleID(iEvent,nCompton)=A(3);
                        Compton.SourceID(iEvent,nCompton)=A(4);
                        Compton.BaseID(iEvent,nCompton)=A(5);
                        Compton.CrystalID(iEvent,nCompton)=A(6);
                        Compton.PixelID(iEvent,nCompton)=A(7);
                        Compton.Time(iEvent,nCompton)=A(8);
                        Compton.Energy(iEvent,nCompton)=A(9);
                        Compton.Range(iEvent,nCompton)=A(10);
                        Compton.XPos(iEvent,nCompton)=A(11);
                        Compton.YPos(iEvent,nCompton)=A(12);
                        Compton.ZPos(iEvent,nCompton)=A(13);
                        Compton.GParticle(iEvent,nCompton)=A(14);
                        Compton.HitParticleID(iEvent,nCompton)=A(15);
                        Compton.MotherParticleID(iEvent,nCompton)=A(16);
                        Compton.PhotonID(iEvent,nCompton)=A(17);
                        Compton.NCompton(iEvent,nCompton)=A(18);
                        Compton.NRayleigh(iEvent,nCompton)=A(19);
                        %
                        %
                        %TotalEnergy=zeros(nEvents,1);
                        %NEventComptons=zeros(nEvents,1);
                        %XFirstHit=zeros(nEvents,2);
                        TotalEnergy(iEvent)=TotalEnergy(iEvent)+(1000.*A(9));
                        NEventComptons(iEvent)=NEventComptons(iEvent)+1;
                        if (nCompton==1)
                          XFirstHit(iEvent)=A(11);
                        end                        
                    case 'phot'
                        nPhot=nPhot+1;
                        Phot.nPhot(iEvent)=nPhot;
                        Phot.nHit(iEvent,nPhot)=nHit;
                        Phot.RunID(iEvent,nPhot)=A(1);
                        Phot.EventID(iEvent,nPhot)=A(2);
                        Phot.PrimaryParticleID(iEvent,nPhot)=A(3);
                        Phot.SourceID(iEvent,nPhot)=A(4);
                        Phot.BaseID(iEvent,nPhot)=A(5);
                        Phot.CrystalID(iEvent,nPhot)=A(6);
                        Phot.PixelID(iEvent,nPhot)=A(7);
                        Phot.Time(iEvent,nPhot)=A(8);
                        Phot.Energy(iEvent,nPhot)=A(9);
                        Phot.Range(iEvent,nPhot)=A(10);
                        Phot.XPos(iEvent,nPhot)=A(11);
                        Phot.YPos(iEvent,nPhot)=A(12);
                        Phot.ZPos(iEvent,nPhot)=A(13);
                        Phot.GParticle(iEvent,nPhot)=A(14);
                        Phot.HitParticleID(iEvent,nPhot)=A(15);
                        Phot.MotherParticleID(iEvent,nPhot)=A(16);
                        Phot.PhotonID(iEvent,nPhot)=A(17);
                        Phot.NCompton(iEvent,nPhot)=A(18);
                        Phot.NRayleigh(iEvent,nPhot)=A(19);
                        %
                        %TotalEnergy=zeros(nEvents,1);
                        %XPhotocapture=zeros(nEvents,1);
                        %YPhotocapture=zeros(nEvents,1);
                        %ZPhotocapture=zeros(nEvents,1);
                        TotalEnergy(iEvent)=TotalEnergy(iEvent)+(1000.*A(9));
                        XPhotocapture(iEvent)=A(11);
                        YPhotocapture(iEvent)=A(12);
                        ZPhotocapture(iEvent)=A(13);
                      case 'Cerenkov'
                        nCerenkov=nCerenkov+1;
                        Cerenkov.nHit(iEvent,nCerenkov)=nHit;
                        Cerenkov.RunID(iEvent,nCerenkov)=A(1);
                        Cerenkov.EventID(iEvent,nCerenkov)=A(2);
                        Cerenkov.PrimaryParticleID(iEvent,nCerenkov)=A(3);
                        Cerenkov.SourceID(iEvent,nCerenkov)=A(4);
                        Cerenkov.BaseID(iEvent,nCerenkov)=A(5);
                        Cerenkov.CrystalID(iEvent,nCerenkov)=A(6);
                        Cerenkov.PixelID(iEvent,nCerenkov)=A(7);
                        Cerenkov.Time(iEvent,nCerenkov)=A(8);
                        Cerenkov.Energy(iEvent,nCerenkov)=A(9);
                        Cerenkov.Range(iEvent,nCerenkov)=A(10);
                        Cerenkov.XPos(iEvent,nCerenkov)=A(11);
                        Cerenkov.YPos(iEvent,nCerenkov)=A(12);
                        Cerenkov.ZPos(iEvent,nCerenkov)=A(13);
                        Cerenkov.GParticle(iEvent,nCerenkov)=A(14);
                        Cerenkov.HitParticleID(iEvent,nCerenkov)=A(15);
                        Cerenkov.MotherParticleID(iEvent,nCerenkov)=A(16);
                        Cerenkov.CerenkovonID(iEvent,nCerenkov)=A(17);
                        Cerenkov.NCompton(iEvent,nCerenkov)=A(18);
                        Cerenkov.NRayleigh(iEvent,nCerenkov)=A(19);
                    case 'msc'
                        nMsc=nMsc+1;
                        Msc.nMsc(iEvent)=nMsc;
                    case 'Rayl'
                        nRayl=nRayl+1;
                        Rayl.nRayl(iEvent)=nRayl;
                    case 'eIoni'
                        neIoni=neIoni+1;
                        eIoni.neIoni(iEvent)=neIoni;
                    case 'eBrem'
                        neBrem=neBrem+1;
                        eBrem.neBrem(iEvent)=neBrem;
                    case 'Scintillation'
                        nScintillation=nScintillation+1;
                        Scintillation.nScintillation(iEvent)=nScintillation;
                    otherwise
                        disp(['Unusual Str1 Content: ',Str1]);
                        disp(tline);
                end
            end
        end
    end
end
fclose(fid);
%
%%%%%%%%%%%%%%%5
%
%
HitEvents.TotalEnergy=TotalEnergy;
HitEvents.NEventComptons=NEventComptons;
HitEvents.XFirstHit=XFirstHit;
HitEvents.XPhotocapture=XPhotocapture;
HitEvents.YPhotocapture=YPhotocapture;
HitEvents.ZPhotocapture=ZPhotocapture;
%
save ReadHits.mat Trans Absorption Compton Phot Cerenkov Msc Rayl eIoni eBrem  ...
  HitEvents Scintillation Total -v7