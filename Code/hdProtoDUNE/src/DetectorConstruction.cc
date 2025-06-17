//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm5/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "Run.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "PMTSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include <map>
#include <string>
#include <iterator>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberMaterial(nullptr),fWorldMaterial(nullptr),
   fSolidWorld(nullptr),fLogicWorld(nullptr),fPhysiWorld(nullptr),
   fSolidAbsorber(nullptr),fLogicAbsorber(nullptr),fPhysiAbsorber(nullptr),
   fDetectorMessenger(nullptr),lvol_fiber(nullptr),pvol_fiber(nullptr)
{
  // default parameter values of the calorimeter
  fAbsorberThickness = 75.*mm;
  fAbsorberSizeYZ    = 47.*mm;
  fXposAbs           = 0.*cm;
  fiberLength = 5.*cm;
  ComputeGeomParameters();
  
  // materials  
  DefineMaterials();
  SetWorldMaterial   ("G4_Galactic");
  SetAbsorberMaterial("G4_lAr");
  //SetAbsorberMaterial("G4_C");

  auto* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
  //auto* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  std::vector<G4double> energy     = {9.69*eV,9.89*eV,10.0*eV};
  // refractive index from Collins [https://indico.cern.ch/event/609440/contributions/2548164/attachments/1441824/2220090/Light_Propagation_in_Liquid_Argon.pdf]
  std::vector<G4double> rindex     = {1.42,1.43,1.44};
  // taking the value of 1ppm C impurities in lAr with cross section 1Mbarn per molecule
  // absorption length falls linearly with impurity concentration , cross section per molecule
  G4double absLen = 1.5*m;
  absorptionLength = {absLen, absLen, absLen};

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  // property independent of energy

  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  tpb = 0.02;
  // Doke et al. , NIMA 291,3(1990) [https://www.sciencedirect.com/science/article/pii/016890029090011T?via%3Dihub] 
  MPT->AddConstProperty("SCINTILLATIONYIELD", tpb*51300./MeV);
  MPT->AddProperty("RINDEX", energy, rindex);
  MPT->AddProperty("ABSLENGTH", energy, absorptionLength);

  // fano factor of LAr = 0.11 , Doke et al, NIM 134 (1976)353, taken from [https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf]

  // The actual number of emitted photons during a step
  //fluctuates around the mean number of photons with
  //a width given by ResolutionScale x sqrt(MeanNumberOfPhotons)
  MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);

  MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 7. * ns);
  // we take an average of the following for the slow time components:
  // [https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf] - 1400ns
  // [https://indico.cern.ch/event/609440/contributions/2548164/attachments/1441824/2220090/Light_Propagation_in_Liquid_Argon.pdf] - 1500ns
  // [Creus et al. , JINST 10 (2015) ] - 1600ns
  MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 1500. * ns);
  // values taken from [https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf] - yield fraction = 0.75 in old geant4 version
  MPT->AddConstProperty("SCINTILLATIONYIELD1", 0.75);
  MPT->AddConstProperty("SCINTILLATIONYIELD2", 0.25);
  MPT->AddProperty("SCINTILLATIONCOMPONENT1", energy, { 1.0, 1.0, 1.0 });
  MPT->AddProperty("SCINTILLATIONCOMPONENT2", energy, { 1.0, 1.0, 1.0 });

  lAr->SetMaterialPropertiesTable(MPT);

  // Set the Birks Constant for the LAr scintillator, using C. William thesis [http://unizh.web.cern.ch/Publications/Articles/thesis_William.pdf]
  //lAr->GetIonisation()->SetBirksConstant(0.053* mm / MeV);
  // trying using Doke et al. 1990 paper
 
 //lAr->GetIonisation()->SetBirksConstant(0.89* mm / MeV);
 birks = 0.03*mm /MeV;
 lAr->GetIonisation()->SetBirksConstant(birks);
  std::cout << "lAr Ion Pair Mean Energy = " << lAr->GetIonisation()->GetMeanEnergyPerIonPair()/CLHEP::eV << " eV" << std::endl;
  if(fLogicAbsorber) { fLogicAbsorber->SetMaterial(lAr); }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

  //SetAbsorberMaterial("G4_lAr");
 //   fLogicAbsorber->SetMaterial(lAr);
    //G4RunManager::GetRunManager()->PhysicsHasBeenModified();
 
  // create commands for interactive definition of the calorimeter  
  fDetectorMessenger = new DetectorMessenger(this,run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  //This function illustrates the possible ways to define materials
 
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  

  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;
  
  //
  // define Elements
  //

  G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);
  G4Element* Ar = new G4Element("Argon",   symbol="Ar", z=18, a=  39.95*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I" , z=53, a= 126.90*g/mole);
  G4Element* Xe = new G4Element("Xenon",   symbol="Xe", z=54, a= 131.29*g/mole);

  Al = new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);

  G4Material* lAr = 
    new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
  lAr->AddElement(Ar, natoms=1);

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  // Scintillator

  G4Material* Sci =
       new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C, natoms=8);
  Sci->AddElement(H, natoms=8);
     
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  std::vector<G4double> energy     = {9.5*eV,9.6*eV,9.7*eV};
  std::vector<G4double> rindex     = {1.42,1.43,1.44};
  std::vector<G4double> absorptionLength2 = {150.*cm, 150.*cm, 150.*cm};

  G4MaterialPropertiesTable* MPT2 = new G4MaterialPropertiesTable();


  MPT2->AddConstProperty("SCINTILLATIONYIELD", 12000./MeV);
  MPT2->AddProperty("RINDEX", energy, rindex);
  MPT2->AddProperty("ABSLENGTH", energy, absorptionLength2);


  MPT2->AddConstProperty("RESOLUTIONSCALE",1.0);

  MPT2->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 7. * ns);
  MPT2->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  MPT2->AddProperty("SCINTILLATIONCOMPONENT1", energy, { 1.0, 1.0, 1.0 });

  Sci->SetMaterialPropertiesTable(MPT2);
  fiberMat = Sci;

  //
  // example of vacuum
  //
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                 kStateGas,temperature,pressure);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeGeomParameters()
{
  // Compute derived parameters of the calorimeter
  fXstartAbs = fXposAbs-0.5*fAbsorberThickness; 
  fXendAbs   = fXposAbs+0.5*fAbsorberThickness;

  G4double xmax = std::max(std::abs(fXstartAbs), std::abs(fXendAbs)) +0.73*m+fiberLength;
  std::cout << "xmaxm = " << xmax << std::endl;
 // fWorldSizeX = 2.4*xmax; 
  //fWorldSizeYZ= fWorldSizeX;
  fWorldSizeX = fWorldSizeYZ =  20.*m;

  if(nullptr != fPhysiWorld) { ChangeGeometry(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  if(nullptr != fPhysiWorld) { return fPhysiWorld; }
  // World
  //
  std::cout << "xxxxxxxxxxxxxxxxxxx" << " world size = " << G4ThreeVector(fWorldSizeX,fWorldSizeYZ,fWorldSizeYZ)/CLHEP::cm << " cm" << std::endl; 

  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);  //its size
                         
  fLogicWorld = new G4LogicalVolume(fSolidWorld,          //its solid
                                   fWorldMaterial,        //its material
                                   "World");              //its name
                                   
  fPhysiWorld = new G4PVPlacement(0,                      //no rotation
                                 G4ThreeVector(0.,0.,0.), //at (0,0,0)
                                 fLogicWorld,             //its logical volume
                                 "World",                 //its name
                                 0,                       //its mother  volume
                                 false,                   //no boolean operation
                                 0);                      //copy number
   
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;
  
  G4int ncomponents, natoms;
  G4double fractionmass;

  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4double airSizeX = fWorldSizeX/2;
  G4double airSizeYZ = fWorldSizeYZ/2;


  airSolid = new G4Box("air",
                      airSizeX,airSizeYZ,airSizeYZ);                              
  airLogic = new G4LogicalVolume(airSolid,    //its solid
                                Air,
                                "air");

  airPhysi = new G4PVPlacement(0,                      //no rotation
                              G4ThreeVector(0,0.,0.), 
                              airLogic,             //its logical volume
                             "air",                 //its name
                              fLogicWorld,                       //its mother  volume
                             false,                   //no boolean operation
                             0);


 
  auto* Teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
  std::vector<G4double> energy     = {9.69*eV,9.89*eV,10.0*eV};
  std::vector<G4double> rindex     = {1.35,1.37,1.38};
  G4MaterialPropertiesTable* MPT3 = new G4MaterialPropertiesTable();

  MPT3->AddProperty("RINDEX", energy, rindex);
  Teflon->SetMaterialPropertiesTable(MPT3);

  std::cout << "duuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuummmmmmmmmmmmmmmmmmmmmmmmmmmmmmpppppppppppppppppppppppppppiiiiiiiiiiiiiiiiiiiiiiiiinnnnng" << std::endl;

  // ARGON CELL
  LArX = 4.6*m;
  LArY = 6.1*m;
  LArZ = 7.2*m;

  std::cout << "after LArX def" << std::endl;

  Teflon->GetMaterialPropertiesTable()->DumpTable();

  //fSolidTeflon = new G4Box("Teflon",0.5*LArX+1.*cm,0.5*LArY+1.*cm,0.5*LArZ+1.*cm);
    
  //fLogicTeflon = new G4LogicalVolume(fSolidTeflon,    //its solid
  //                             Teflon,                //its material
    //                             "Teflon");       //its name
    
  G4ThreeVector u = G4ThreeVector(1 , 0, 0);         
  G4ThreeVector v = G4ThreeVector(0 , 1, 0);        
  G4ThreeVector w = G4ThreeVector(0, 0, 1);
    
  G4RotationMatrix rotm1  = G4RotationMatrix(u, v, w);
  G4ThreeVector position1 = G4ThreeVector(0,0,0);
  const G4Transform3D & transform1 = G4Transform3D(rotm1,position1);
    
  //fPhysiTeflon = new G4PVPlacement(transform1,
    //                 fLogicTeflon,     //its logical volume
      //                 "Teflon",         //its name
        //                 airLogic,        //its mother
          //                 false,              //no boulean operat
            //                 0);                 //copy number

  fSolidAbsorber = new G4Box("lAr",0.5*LArX,0.5*LArY,0.5*LArZ);
                          
  fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber,    //its solid
               fAbsorberMaterial, //its material
                 "Absorber");       //its name
 
  fPhysiAbsorber = new G4PVPlacement(0,G4ThreeVector(),
      fLogicAbsorber,     //its logical volume
        "Absorber",         //its name
          airLogic,        //its mother
            false,              //no boulean operat
              0);                 //copy number
 

  /////////////////////////////////////////////////////////////////    SENSITIVE PMTS ON LAr BASES ///////////////////////////////////////////////////////////////////////////////////////////
  // Sensitive Detectors on Cylinder Bases (representing PMTs)
  //G4Material* pmtMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  //G4Tubs* pmtSurface = new G4Tubs("PMTSurface", 0., lArDiameter/2, 0.5 * mm, 0., 360. * deg);

  // Top Base
  //G4ThreeVector positionTop = G4ThreeVector(0,0,lArLength/2 - 1.*mm);
  //const G4Transform3D & transformTop = G4Transform3D(rotm1,positionTop);
  //G4LogicalVolume* pmtTopLog = new G4LogicalVolume(pmtSurface, pmtMaterial, "PMT_Top");
  //new G4PVPlacement(nullptr,positionTop, pmtTopLog, "PMT_Top", fLogicAbsorber, false, 0);

  // Bottom Base
  //G4ThreeVector positionBottom = G4ThreeVector(0,0,-lArLength/2+1.*mm);
  //const G4Transform3D & transformBottom = G4Transform3D(rotm1,positionBottom);
  //G4LogicalVolume* pmtBottomLog = new G4LogicalVolume(pmtSurface, pmtMaterial, "PMT_Bottom");
  //new G4PVPlacement(nullptr,positionBottom, pmtBottomLog, "PMT_Bottom", fLogicAbsorber, false, 0);

  // Attach Sensitive Detector to Bases
  //G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  //EventAction* eventAction = dynamic_cast<EventAction*>(
    //const_cast<G4UserEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction()));


  //auto pmtSD = new PMTSensitiveDetector("PMT_SD",eventAction);
  //sdManager->AddNewDetector(pmtSD);
  //pmtTopLog->SetSensitiveDetector(pmtSD);
  //pmtBottomLog->SetSensitiveDetector(pmtSD);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //G4OpticalSurface* OpSurface = new G4OpticalSurface("name");
  //G4LogicalBorderSurface* Surface = new G4LogicalBorderSurface("name",fPhysiAbsorber,fPhysiTeflon,OpSurface);

  //OpSurface->SetType(dielectric_metal);
  //OpSurface->SetModel(unified);
  //OpSurface->SetFinish(polished);
  //OpSurface->SetSigmaAlpha(0.1);

  //std::vector<G4double> pp = {9.6*eV, 10.0*eV};
  //std::vector<G4double> specularlobe = {0.3, 0.3};
  //std::vector<G4double> specularspike = {0.2, 0.2};
  //std::vector<G4double> backscatter = {0.1, 0.1};
  //rindex = {1.35, 1.40};
  //std::vector<G4double> reflectivity = {0.98, 0.98};
  //std::vector<G4double> efficiency = {1.0, 1.0};

  //G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();

  //SMPT->AddProperty("RINDEX", pp, rindex);
  //SMPT->AddProperty("SPECULARLOBECONSTANT", pp, specularlobe);
  //SMPT->AddProperty("SPECULARSPIKECONSTANT", pp, specularspike);
  //SMPT->AddProperty("BACKSCATTERCONSTANT", pp, backscatter);
  //SMPT->AddProperty("REFLECTIVITY", pp, reflectivity);
  //SMPT->AddProperty("EFFICIENCY", pp, efficiency);

  //OpSurface->SetMaterialPropertiesTable(SMPT);
 

  // ORGANIC SCINTILLATOR
  //std::cout << "fiber length = " << fiberLength/CLHEP::cm << " cm" << std::endl;
  //std::cout << "fiber diameter = " << fiberDiameter/CLHEP::cm << " cm" << std::endl;
  
  //svol_fiber = new G4Tubs("fiber",                      //name
                         //0*mm, 0.5*fiberDiameter,       //r1, r2
                         //0.5*fiberLength,               //half-length
                         //0., twopi);                    //theta1, theta2

  //lvol_fiber = new G4LogicalVolume(svol_fiber,          //solid
    //fiberMat,            //material
      //"fiber");            //name


  // rotation around y axis
  // u, v, w are the daughter axes, projected on the mother frame
  //G4double rotationAngle = pi/2 - pi/2;
  //u = G4ThreeVector(std::cos(rotationAngle) , 0, std::sin(rotationAngle));
  //v = G4ThreeVector(0 , 1, 0);
  //w = G4ThreeVector(std::sin(-rotationAngle), 0, std::cos(rotationAngle)); 


  // rotation around x axis
  //u = G4ThreeVector(1,0,0);
  //v = G4ThreeVector(0,std::cos(rotationAngle),std::sin(-rotationAngle));
  //w = G4ThreeVector(0,std::sin(rotationAngle),std::cos(rotationAngle));
 
  //G4RotationMatrix rotm2  = G4RotationMatrix(u, v, w);
  //G4ThreeVector position2 = G4ThreeVector(1.3*m,0.*cm,0);

  //std::cout << "position before rotation = " << position2 << std::endl;
  //position2 = rotateAroundY(position2,pi/2);
  //std::cout << "position after rotation = " << position2 << std::endl;
  //const G4Transform3D & transform2 = G4Transform3D(rotm2,position2);


 //pvol_fiber = new G4PVPlacement(transform2,
                                //lvol_fiber,                    //its logical volume
                                //"liquidScintillator",          //its name
                                //airLogic,                      //its mother
                                //false,                         //no boulean operat
                                //0); 

 std::vector<G4VisAttributes*> fVisAttributes;
 auto visAttributes = new G4VisAttributes(G4Colour(0.,0.,0.)); // world is black
 //visAttributes->SetVisibility(false);
 fLogicWorld->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);                                     

 visAttributes = new G4VisAttributes(G4Colour(0.4,0.4,0.4));   // make inner LAr detector dark grey
 fLogicAbsorber->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,0.0));   // make Aluminium Layer yellow
 //fLogicTeflon->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(0.4,0.4,0.4,0.2));   // make outer LAr detector dark grey and opaque
 //biggerLogicAbsorber->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(0.19,0.83,0.78));   // make liquid scintilator turquoise
 //new_lvol_fiber->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(0.19,0.83,0.78));
 //testLogic->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 // PrintGeomParameters();         
  
  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintGeomParameters()
{
  G4cout << "\n" << fWorldMaterial    << G4endl;
  G4cout << "\n" << fAbsorberMaterial << G4endl;
    
  G4cout << "\n The  WORLD   is made of "  << G4BestUnit(fWorldSizeX,"Length")
         << " of " << fWorldMaterial->GetName();
  G4cout << ". The transverse size (YZ) of the world is " 
         << G4BestUnit(fWorldSizeYZ,"Length") << G4endl;
  G4cout << " The ABSORBER is made of " 
         <<G4BestUnit(fAbsorberThickness,"Length")
         << " of " << fAbsorberMaterial->GetName();
  G4cout << ". The transverse size (YZ) is " 
         << G4BestUnit(fAbsorberSizeYZ,"Length") << G4endl;
  G4cout << " X position of the middle of the absorber "
         << G4BestUnit(fXposAbs,"Length");
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  std::cout << "argon40 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,40) << std::endl;
  std::cout << "argon36 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,36) << std::endl;
  std::cout << "argon38 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,38) << std::endl;
  std::cout << "argon39 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,39) << std::endl;
  std::cout << "argon41 abundance = " << G4NistManager::Instance()->GetIsotopeAbundance(18,41) << std::endl;

  if (pttoMaterial && fAbsorberMaterial != pttoMaterial) {
    fAbsorberMaterial = pttoMaterial;

    if(fLogicAbsorber) { fLogicAbsorber->SetMaterial(fAbsorberMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial; 
    if(fLogicWorld) { fLogicWorld->SetMaterial(fWorldMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  fAbsorberThickness = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberSizeYZ(G4double val)
{
  fAbsorberSizeYZ = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeX(G4double val)
{
  fWorldSizeX = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeYZ(G4double val)
{
  fWorldSizeYZ = val;
  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberXpos(G4double val)
{
  fXposAbs = val;
  ComputeGeomParameters();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField()
{
  if ( fFieldMessenger.Get() == 0 ) {
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    G4GlobalMagFieldMessenger* msg =
      new G4GlobalMagFieldMessenger(fieldValue);
    //msg->SetVerboseLevel(1);
    G4AutoDelete::Register(msg);
    fFieldMessenger.Put( msg );        
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ChangeGeometry()
{
  fSolidWorld->SetXHalfLength(fWorldSizeX*0.5);
  fSolidWorld->SetYHalfLength(fWorldSizeYZ*0.5);
  fSolidWorld->SetZHalfLength(fWorldSizeYZ*0.5);

 // fSolidAbsorber->SetXHalfLength(fAbsorberThickness*0.5);
  //fSolidAbsorber->SetYHalfLength(fAbsorberSizeYZ*0.5);
  //fSolidAbsorber->SetZHalfLength(fAbsorberSizeYZ*0.5);

  svol_fiber->SetOuterRadius(fiberDiameter*0.5);
  svol_fiber->SetZHalfLength(fiberLength*0.5);
  svol_fiber->SetDeltaPhiAngle(2*M_PI);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

