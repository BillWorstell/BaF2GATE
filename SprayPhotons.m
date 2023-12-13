% SprayPhotons
%
PhotonsPerMeV=900;
Pi=acos(-1.);
XMin=-25.;
XMax=25.;
%
load ReadHits.mat;
%
NEvents=numel(Compton.nCompton);
%
HiX.YMean=zeros(NEvents,1);
HiX.YStd=zeros(NEvents,1);
HiX.ZMean=zeros(NEvents,1);
HiX.ZStd=zeros(NEvents,1);
LoX.YMean=zeros(NEvents,1);
LoX.YStd=zeros(NEvents,1);
LoX.ZMean=zeros(NEvents,1);
LoX.ZStd=zeros(NEvents,1);
%
for iEvent=1:NEvents
  if (mod(iEvent,100)==0) 
    disp(['iEvent= ',num2str(iEvent)]);
  end
  %
  ThisNCompton=Compton.nCompton(iEvent);
  if (ThisNCompton==0)
    FirstX(iEvent)=Phot.XPos(iEvent,1);
  end
  %
  % Photons from Photocapture interaction
  NPhotonsPhotocapture=round(PhotonsPerMeV*Phot.Energy(iEvent,1));
  NPhotonsComptons=sum(round(PhotonsPerMeV*Compton.Energy(iEvent,1:ThisNCompton)));
  ThisNPhotons=NPhotonsPhotocapture+NPhotonsComptons;
  %
  % Set up array to hold photon detection coordinates
  PhotonX=-zeros(ThisNPhotons,1);
  PhotonY=zeros(ThisNPhotons,1);
  PhotonZ=zeros(ThisNPhotons,1);
  %
  % Position of photocapture interaction
  X0=Phot.XPos(iEvent,1);
  Y0=0.;
  Z0=0.;
  %
  % Isotropically emit scintillation photons from photocapture
  Phi=2.*Pi*rand(NPhotonsPhotocapture,1);
  SinTheta=-1.+2.*rand(NPhotonsPhotocapture,1);
  CosTheta=sqrt(1.-(SinTheta.*SinTheta));
  vX=SinTheta;
  vY=cos(Phi).*CosTheta;
  vZ=sin(Phi).*CosTheta;
  % Forward scintillation
  fvXPos=find(vX>0);
  PhotonT=(XMax-X0)*(1./vX(fvXPos));
  PhotonX(fvXPos)=X0+vX(fvXPos).*PhotonT;
  PhotonY(fvXPos)=Y0+vY(fvXPos).*PhotonT;
  PhotonZ(fvXPos)=Z0+vZ(fvXPos).*PhotonT;
  % Backward scintillation
  fvXNeg=find(vX<0);
  PhotonT=(XMin-X0)*(1./vX(fvXNeg));
  PhotonX(fvXNeg)=X0+vX(fvXNeg).*PhotonT;
  PhotonY(fvXNeg)=Y0+vY(fvXNeg).*PhotonT;
  PhotonZ(fvXNeg)=Z0+vZ(fvXNeg).*PhotonT;
  %
  NPhotonsSoFar=NPhotonsPhotocapture;
  %
  if (ThisNCompton>0)
    
  % Now loop over scintillation photons from Comptons
  for ThisCompton=1:ThisNCompton
    ThisNPhotonsCompton=round(PhotonsPerMeV*Compton.Energy(iEvent,ThisCompton));
    %
    % Position of Compton interaction
    X0=Compton.XPos(iEvent,1);
    Y0=Compton.XPos(iEvent,2);
    Z0=Compton.XPos(iEvent,3);
    %
    % Isotropically emit scintillation photons from Compton interaction
    Phi=2.*Pi*rand(ThisNPhotonsCompton,1);
    SinTheta=-1.+2.*rand(ThisNPhotonsCompton,1);
    CosTheta=sqrt(1.-(SinTheta.*SinTheta));
    vX=SinTheta;
    vY=cos(Phi).*CosTheta;
    vZ=sin(Phi).*CosTheta;
    % Forward scintillation
    fvXPos=find(vX>0);
    PhotonT=(XMax-X0)*(1./vX(fvXPos));
    PhotonX(NPhotonsSoFar+fvXPos)=X0+vX(fvXPos).*PhotonT;
    PhotonY(NPhotonsSoFar+fvXPos)=Y0+vY(fvXPos).*PhotonT;
    PhotonZ(NPhotonsSoFar+fvXPos)=Z0+vZ(fvXPos).*PhotonT;
    % Backward scintillation
    fvXNeg=find(vX<0);
    PhotonT=(XMin-X0)*(1./vX(fvXNeg));
    PhotonX(NPhotonsSoFar+fvXNeg)=X0+vX(fvXNeg).*PhotonT;
    PhotonY(NPhotonsSoFar+fvXNeg)=Y0+vY(fvXNeg).*PhotonT;
    PhotonZ(NPhotonsSoFar+fvXNeg)=Z0+vZ(fvXNeg).*PhotonT;   
    % 
    NPhotonsSoFar=NPhotonsSoFar+ThisNPhotonsCompton;
  end
  end
  %
  % Store Mean and RMS for Y and Z distributions at Low and High X
  fHiX=find(PhotonX>=XMax);
  if (numel(fHiX)>0)
    HiX.YMean(iEvent)=mean(PhotonY(fHiX));
    HiX.YStd(iEvent)=std(PhotonY(fHiX));
    HiX.ZMean(iEvent)=mean(PhotonZ(fHiX));
    HiX.ZStd(iEvent)=std(PhotonZ(fHiX));
  end
  %
  fLoX=find(PhotonX<=XMin);
  if (numel(fLoX)>0)
    LoX.YMean(iEvent)=mean(PhotonY(fLoX));
    LoX.YStd(iEvent)=std(PhotonY(fLoX));
    LoX.ZMean(iEvent)=mean(PhotonZ(fLoX));
    LoX.ZStd(iEvent)=std(PhotonZ(fLoX));
  end
  %
  disp('here');
  %
end

