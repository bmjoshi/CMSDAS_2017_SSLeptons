
class Lepton {
public:
  
  //Four vector
  double pt;
  double eta;
  double phi;

  int charge;

  bool isEl;
  bool isMu;

  bool Tight;
  bool Loose;
};


class Electron : public Lepton{
public:
  
  //variables for tracking cuts
  double dEta; 
  double dPhi;
  //variables for primary vtx cuts
  double dZ;
  double d0;
  // Hadronic energy divided by the Electromagnetic
  double hOverE;
  //Minimum hits
  double mHits;
  double ooEmooP;
  //Isolation
  double relIso;
  double sigmaIetaIeta;
  int    chargeConsistency;

  bool tight(){
    //Barrel
    if(fabs(eta) <= 1.479){
      if(fabs(dEta) > 0.0091)    return false;
      if(fabs(dPhi) > 0.031)     return false;
      if(sigmaIetaIeta > 0.0106) return false;
      if(hOverE > 0.0532)        return false;
      if(fabs(d0) > 0.0126)      return false;
      if(fabs(dZ) > 0.0116)      return false;
      if(fabs(ooEmooP) > 0.0609) return false;
      if(relIso > 0.1649)        return false;
      if(mHits > 1)              return false;
      if(chargeConsistency < 1)  return false;
      if(pt <20)                 return false;
    }
    //Endcap
    else{
      if(fabs(dEta) > 0.0106)    return false;
      if(fabs(dPhi) > 0.0359)    return false;
      if(sigmaIetaIeta > 0.0305) return false;
      if(hOverE > 0.0835)        return false;
      if(fabs(d0) > 0.0163)      return false;
      if(fabs(dZ) > 0.5999)      return false;
      if(fabs(ooEmooP) > 0.1126) return false;
      if(relIso > 0.2075)        return false;
      if(mHits > 1)              return false;
      if(chargeConsistency < 1)  return false;
      if(pt <20)                 return false;
    }
    return true;
  }
  bool loose(){
    //Barrel
    if(fabs(eta) <= 1.479){
      if(fabs(dEta) > 0.0181)    return false;
      if(fabs(dPhi) > 0.0936)    return false;
      if(sigmaIetaIeta > 0.0123) return false;
      if(hOverE > 0.141)         return false;
      if(fabs(d0) > 0.0166)      return false;
      if(fabs(dZ) > 0.54342)     return false;
      if(fabs(ooEmooP) > 0.1353) return false;
      if(relIso > 0.24)          return false;
      if(mHits > 1)              return false;
      if(chargeConsistency < 1)  return false;
      if(pt <20)                 return false;
    }

    //Endcap
    else{
      if(fabs(dEta) > 0.0124)    return false;
      if(fabs(dPhi) > 0.0642)    return false;
      if(sigmaIetaIeta > 0.035)  return false;
      if(hOverE > 0.1115)        return false;
      if(fabs(d0) > 0.098)       return false;
      if(fabs(dZ) > 0.9187)      return false;
      if(fabs(ooEmooP) > 0.1443) return false;
      if(relIso > 0.3529)        return false;
      if(mHits > 1)              return false;
      if(chargeConsistency < 1)  return false;
      if(pt <20)                 return false;
    }
    return true;
  }

  void init(){
    pt   = -100;
    eta  = -100;
    phi  = -100;
    
    dEta = -100;
    dPhi = -100;

    dZ   = -100;
    d0   = -100;
    
    hOverE  = -100;
    
    mHits   =  100;
    ooEmooP = -100;
    
    relIso = 100;
    sigmaIetaIeta = -100;
    
    chargeConsistency = 0;

  }

};

class Muon : public Lepton{
public:

  int isLoose;
  int isTight;

  int    global;
  double chi2;
  int    nValMuHits;
  int    nMatchedStations;
  double dxy;
  double dz;
  int    nValPixelHits;
  int    nTrackerLayers;
  double relIso;


  bool tight(){
    if(!isTight) return false;
    if(pt<20) return false;
    else return true;
  }

  bool loose(){
    if(!isLoose) return false;
    if(pt<20) return false;
    else return true;
  }

  bool cutBasedTight(){
    if (pt < 20)              return false;
    if (global < 5)           return false; //Global muon and tracker muon
    if (chi2 > 10)            return false;
    if (nValMuHits < 1)       return false;
    if (nMatchedStations < 2) return false;
    if (dxy > 0.2)            return false;
    if (dz > 0.5)             return false;
    if (nValPixelHits < 1)    return false;
    if (nTrackerLayers < 6)   return false;
    if (relIso > 0.06)        return false;

    return true;
  }

  bool cutBasedLoose(){
    if (pt < 20)              return false;
    //if (global < 5)           return false; //Global muon and tracker muon
 // if (chi2 > 50)            return false; //Cut relaxed from 10 to 50
    //if (nValMuHits < 1)       return false;
 // if (nMatchedStations < 2) return false; //Cut removed
    if (dxy > 2.0)            return false; //Cut relaxed from 0.2 to 2.0
 // if (dz > 0.5)             return false; //Cut removed 
    if (nValPixelHits < 1)    return false;
 // if (nTrackerLayers < 6)   return false; //Cut removed
    //if (relIso > 0.20)        return false; //Cut relaxed from 0.12 to 0.40

    return true;
  }

  void init(){
    pt   = -100;
    eta  = -100;
    phi  = -100;

    charge = 10.5;

    isLoose = -1;
    isTight = -1;

    global           = -100;
    chi2             =  100;
    nValMuHits       = -100;
    nMatchedStations = -100;
    dxy              =  100;
    dz               =  100;
    nValPixelHits    = -100;
    nTrackerLayers   = -100;
    relIso           =  100;
  }
};
