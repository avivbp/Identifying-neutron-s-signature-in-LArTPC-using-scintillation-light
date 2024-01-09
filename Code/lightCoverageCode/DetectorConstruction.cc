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

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
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
  fiberDiameter = 12.7*cm;
  fiberLength = fiberDiameter;
  ComputeGeomParameters();
  
  // materials  
  DefineMaterials();
  SetWorldMaterial   ("G4_Galactic");
  SetAbsorberMaterial("G4_lAr");
  //SetAbsorberMaterial("G4_C");

  auto* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
  //auto* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  std::vector<G4double> energy     = {9.691*eV,9.692*eV,9.693*eV};
  // refractive index from Collins [https://indico.cern.ch/event/609440/contributions/2548164/attachments/1441824/2220090/Light_Propagation_in_Liquid_Argon.pdf]
  std::vector<G4double> rindex     = {1.42,1.43,1.44};
  // taking the value of 1ppm C impurities in lAr with cross section 1Mbarn per molecule
  // absorption length falls linearly with impurity concentration , cross section per molecule
  std::vector<G4double> absorptionLength = {50.*cm, 50.*cm, 50.*cm};

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  // property independent of energy

  // Doke et al. , NIMA 291,3(1990) [https://www.sciencedirect.com/science/article/pii/016890029090011T?via%3Dihub] 
  MPT->AddConstProperty("SCINTILLATIONYIELD", 51300./MeV);
  MPT->AddProperty("RINDEX", energy, rindex);
  MPT->AddProperty("ABSLENGTH", energy, absorptionLength);

  // fano factor of LAr = 0.11 , Doke et al, NIM 134 (1976)353, taken from [https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf]

  // The actual number of emitted photons during a step
  //fluctuates around the mean number of photons with
  //a width given by ResolutionScale x sqrt(MeanNumberOfPhotons)
  MPT->AddConstProperty("RESOLUTIONSCALE", 0.11);

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
  lAr->GetIonisation()->SetBirksConstant(0.89* mm / MeV);
  std::cout << "lAr Ion Pair Mean Energy = " << lAr->GetIonisation()->GetMeanEnergyPerIonPair()/CLHEP::eV << " eV" << std::endl;
  if(fLogicAbsorber) { fLogicAbsorber->SetMaterial(lAr); }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

  //SetAbsorberMaterial("G4_lAr");
 //   fLogicAbsorber->SetMaterial(lAr);
    //G4RunManager::GetRunManager()->PhysicsHasBeenModified();
 
  // create commands for interactive definition of the calorimeter  
  fDetectorMessenger = new DetectorMessenger(this);
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

  //
  // define simple materials
  //

  // CREATING DEUTERIUM
  G4Isotope* D  = new G4Isotope("Deuteron", 1, 2, 2.0141018* CLHEP::g / CLHEP::mole);
  G4Element* elD = new G4Element("Deuterium","elD", 1);
  elD->AddIsotope(D, 1);
  G4Material* matD = new G4Material("matD", 0.00018* CLHEP::g / CLHEP::cm3, 1);
  matD->AddElement(elD, 1);

  new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
  Al = new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);

  G4Material* lAr = 
    new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
  lAr->AddElement(Ar, natoms=1);

  new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
  new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
  new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
  new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
  new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);

  //
  // define a material from elements.   case 1: chemical molecule
  //

  G4Material* H2O = new G4Material("Water",density= 1.000*g/cm3,ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78*eV);

  G4Material* CH = new G4Material("Plastic",density= 1.04*g/cm3,ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);

  G4Material* NaI = new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Air20 = 
    new G4Material("Air20", density= 1.205*mg/cm3, ncomponents=2,
                   kStateGas, 293.*kelvin, 1.*atmosphere);
  Air20->AddElement(N, fractionmass=0.7);
  Air20->AddElement(O, fractionmass=0.3);

  //Graphite
  //
  G4Material* Graphite = 
    new G4Material("Graphite", density= 1.7*g/cm3, ncomponents=1);
  Graphite->AddElement(C, fractionmass=1.);

  //Havar
  //
  G4Element* Cr = new G4Element("Chrome", "Cr", z=24, a=  51.996*g/mole);
  G4Element* Fe = new G4Element("Iron"  , "Fe", z=26, a=  55.845*g/mole);
  G4Element* Co = new G4Element("Cobalt", "Co", z=27, a=  58.933*g/mole);
  G4Element* Ni = new G4Element("Nickel", "Ni", z=28, a=  58.693*g/mole);
  G4Element* W  = new G4Element("Tungsten","W", z=74, a= 183.850*g/mole);

  G4Material* Havar = 
    new G4Material("Havar", density= 8.3*g/cm3, ncomponents=5);
  Havar->AddElement(Cr, fractionmass=0.1785);
  Havar->AddElement(Fe, fractionmass=0.1822);
  Havar->AddElement(Co, fractionmass=0.4452);
  Havar->AddElement(Ni, fractionmass=0.1310);
  Havar->AddElement(W , fractionmass=0.0631);

  //
  // examples of gas
  //  
  new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 273.15*kelvin, 1*atmosphere);
                           
  new G4Material("XenonGas", z=54, a=131.29*g/mole, density= 5.458*mg/cm3,
                 kStateGas, 293.15*kelvin, 1*atmosphere);
                           
  G4Material* CO2 =
    new G4Material("CarbonicGas", density= 1.977*mg/cm3, ncomponents=2);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  G4Material* ArCO2 =
    new G4Material("ArgonCO2",   density= 1.8223*mg/cm3, ncomponents=2);
  ArCO2->AddElement (Ar,  fractionmass=0.7844);
  ArCO2->AddMaterial(CO2, fractionmass=0.2156);

  //another way to define mixture of gas per volume
  G4Material* NewArCO2 =
    new G4Material("NewArgonCO2", density= 1.8223*mg/cm3, ncomponents=3);
  NewArCO2->AddElement (Ar, natoms=8);
  NewArCO2->AddElement (C,  natoms=2);
  NewArCO2->AddElement (O,  natoms=4);

  G4Material* ArCH4 = 
    new G4Material("ArgonCH4",    density= 1.709*mg/cm3,  ncomponents=3);
  ArCH4->AddElement (Ar, natoms=93);
  ArCH4->AddElement (C,  natoms=7);
  ArCH4->AddElement (H,  natoms=28);

  G4Material* XeCH = 
    new G4Material("XenonMethanePropane", density= 4.9196*mg/cm3, ncomponents=3,
                   kStateGas, 293.15*kelvin, 1*atmosphere);
  XeCH->AddElement (Xe, natoms=875);
  XeCH->AddElement (C,  natoms=225);
  XeCH->AddElement (H,  natoms=700);

  G4Material* steam = 
    new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1);
  steam->AddMaterial(H2O, fractionmass=1.);
  steam->GetIonisation()->SetMeanExcitationEnergy(71.6*eV);  

  G4Material* rock1 = new G4Material("StandardRock",
                                     2.65*CLHEP::g/CLHEP::cm3, 1, kStateSolid);
  rock1->AddElement(Na, 1);

  // Scintillator

  G4Material* Sci =
       new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C, natoms=8);
  Sci->AddElement(H, natoms=8);
     
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
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

  G4double airSizeX = fWorldSizeX/2 - 0.3*m;
  G4double airSizeYZ = fWorldSizeYZ/2 -30.*cm;


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

  // Absorber
  // 

 // fSolidAbsorber = new G4Box("Absorber",        
                     // fAbsorberThickness/2,fAbsorberSizeYZ/2,fAbsorberSizeYZ/2);

  G4double biggerlArDiameter = 141.*mm;
  G4double biggerlArLength = 150.*mm;
 
  biggerSolidAbsorber = new G4Box("lAr",0.5*m,1.95*m,1.94*m);
  
  biggerLogicAbsorber = new G4LogicalVolume(biggerSolidAbsorber,    //its solid
                                       fAbsorberMaterial, //its material
                                       "biggerAbsorber");       //its name
  
  G4ThreeVector u = G4ThreeVector(1 , 0, 0);
  G4ThreeVector v = G4ThreeVector(0 , 0, -1);
  G4ThreeVector w = G4ThreeVector(0, 1, 0);
  
  G4RotationMatrix rotm0  = G4RotationMatrix(u, v, w);
  G4ThreeVector position0 = G4ThreeVector(0.,0.,0.);
  const G4Transform3D & transform0 = G4Transform3D(rotm0,position0);
  
  biggerPhysiAbsorber = new G4PVPlacement(transform0,
                                     biggerLogicAbsorber,     //its logical volume
                                     "biggerAbsorber",         //its name
                                     airLogic,        //its mother
                                     false,              //no boulean operat
                                     0);

  G4double lArDiameter = 75.*mm;
  G4double lArLength = 47.*mm;

  //auto* testSolid = new G4Box("test",1.5*mm,10.*cm,10.*cm);
  //auto* testLogic = new G4LogicalVolume(testSolid,Al,"logicAl");
  //auto* testPhysi = new G4PVPlacement(0,G4ThreeVector(-0.5*m,0,0),testLogic,"PhysAl",airLogic,false,0);


  //fSolidAluminium = new G4Tubs("Aluminium",0.5*lArDiameter,0.5*lArDiameter+1.5*mm,0.5*lArLength,0.,twopi);
    
  //fLogicAluminium = new G4LogicalVolume(fSolidAluminium,    //its solid
                               //         Al,                //its material
                                 //      "Aluminium");       //its name
    
    //u = G4ThreeVector(1 , 0, 0);         
    //v = G4ThreeVector(0 , 0, -1);        
    //w = G4ThreeVector(0, 1, 0);
    
    //G4RotationMatrix rotm1  = G4RotationMatrix(u, v, w);
    //G4ThreeVector position1 = G4ThreeVector(0,0,0);
    //const G4Transform3D & transform1 = G4Transform3D(rotm1,position1);
    
    //fPhysiAluminium = new G4PVPlacement(0,
                   //                    G4ThreeVector(0,0,0),
                     //                  fLogicAluminium,     //its logical volume
                       //                "Aluminium",         //its name
                         //              biggerLogicAbsorber,        //its mother
                           //            false,              //no boulean operat
                             //          0);                 //copy number


  //fSolidAbsorber = new G4Tubs("lAr",0*cm,0.5*lArDiameter,0.5*lArLength,0.,twopi);
                          
  //fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber,    //its solid
               //                        fAbsorberMaterial, //its material
                 //                      "Absorber");       //its name
 
  //u = G4ThreeVector(1 , 0, 0);
  //v = G4ThreeVector(0 , 0, -1);
  //w = G4ThreeVector(0, 1, 0);
   
  //G4RotationMatrix rotm1  = G4RotationMatrix(u, v, w);
  //G4ThreeVector position1 = G4ThreeVector(0,0,0);
  //const G4Transform3D & transform1 = G4Transform3D(rotm1,position1);
                                            
  //fPhysiAbsorber = new G4PVPlacement(0,
    //                            G4ThreeVector(0,0,0),
      //                          fLogicAbsorber,     //its logical volume
        //                        "Absorber",         //its name
          //                      biggerLogicAbsorber,        //its mother
            //                    false,              //no boulean operat
              //                  0);                 //copy number
 

 // G4OpticalSurface* OpSurface = new G4OpticalSurface("name");

  //G4LogicalBorderSurface* Surface = new
    //G4LogicalBorderSurface("name",fPhysiAbsorber,fPhysiAluminium,OpSurface);

  //OpSurface->SetType(dielectric_dielectric);
  //OpSurface->SetModel(unified);
  //OpSurface->SetFinish(groundfrontpainted);
  //OpSurface->SetSigmaAlpha(0.1);

  //std::vector<G4double> pp = {9.6*eV, 10.0*eV};
  //std::vector<G4double> specularlobe = {0.3, 0.3};
  //std::vector<G4double> specularspike = {0.2, 0.2};
  //std::vector<G4double> backscatter = {0.1, 0.1};
  //std::vector<G4double> rindex = {1.35, 1.40};
  //std::vector<G4double> reflectivity = {0.97, 0.98};
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
  std::cout << "fiber length = " << fiberLength/CLHEP::cm << " cm" << std::endl;
  std::cout << "fiber diameter = " << fiberDiameter/CLHEP::cm << " cm" << std::endl;
  
  svol_fiber = new G4Tubs("fiber",                      //name
                         0*mm, 0.5*fiberDiameter,       //r1, r2
                         0.5*fiberLength,               //half-length
                         0., twopi);                    //theta1, theta2

  lvol_fiber = new G4LogicalVolume(svol_fiber,          //solid
    fiberMat,            //material
      "fiber");            //name


  // rotation around y axis
  // u, v, w are the daughter axes, projected on the mother frame
  G4double rotationAngle = pi/2 - 0*pi/2;
  //u = G4ThreeVector(std::cos(rotationAngle) , 0, std::sin(rotationAngle));
  //v = G4ThreeVector(0 , 1, 0);
  //w = G4ThreeVector(std::sin(-rotationAngle), 0, std::cos(rotationAngle)); 


  // rotation around x axis
  u = G4ThreeVector(1,0,0);
  v = G4ThreeVector(0,std::cos(rotationAngle),std::sin(-rotationAngle));
  w = G4ThreeVector(0,std::sin(rotationAngle),std::cos(rotationAngle));
 
  G4RotationMatrix rotm2  = G4RotationMatrix(u, v, w);
  // rotm2 = G4RotationMatrix(G4ThreeVector(1,0,0),G4ThreeVector(0,1,0),G4ThreeVector(0,0,1));
  G4ThreeVector position2 = G4ThreeVector(-30.17*cm,-2.0*m,0);
  //position2 = G4ThreeVector(20.17*cm,-2.6118*m,0);
  position2 = G4ThreeVector(-58.7*cm,-200.3*cm,0);
  std::cout << "position before rotation = " << position2 << std::endl;
  //position2 = rotateAroundY(position2,0*pi/3);
  std::cout << "position after rotation = " << position2 << std::endl;
  const G4Transform3D & transform2 = G4Transform3D(rotm2,position2);


 pvol_fiber = new G4PVPlacement(transform2,                   
                                lvol_fiber,                    //its logical volume
                                "liquidScintillator",          //its name
                                airLogic,                      //its mother
                                false,                         //no boulean operat
                                0); 

  //new_svol_fiber = new G4Tubs("fiber2",                      //name
                              //0*mm, 0.5*fiberDiameter,       //r1, r2
                              //0.5*fiberLength,               //half-length
                              //0., twopi);                    //theta1, theta2
    
  //new_lvol_fiber = new G4LogicalVolume(new_svol_fiber,          //solid
                                    //fiberMat,            //material
                                    //"fiber2");            //name
    
    
  // u, v, w are the daughter axes, projected on the mother frame
  rotationAngle = pi/2;
  u = G4ThreeVector(std::cos(rotationAngle) , 0, std::sin(rotationAngle));
  v = G4ThreeVector(0 , 1, 0);
  w = G4ThreeVector(std::sin(-rotationAngle), 0, std::cos(rotationAngle));
    
  rotm2  = G4RotationMatrix(u, v, w);
  // rotm2 = G4RotationMatrix(G4ThreeVector(1,0,0),G4ThreeVector(0,1,0),G4ThreeVector(0,0,1));
  //position2 = G4ThreeVector(-0.5*m,0,0);
  std::cout << "position before rotation = " << position2 << std::endl;
  position2 = rotateAroundY(position2,0);
  std::cout << "position after rotation = " << position2 << std::endl;
  const G4Transform3D & transform3 = G4Transform3D(rotm2,position2);
    
    
  //new_pvol_fiber = new G4PVPlacement(transform3,
                                 //new_lvol_fiber,                     //its logical volume
                                 //"liquidScintillator2",           //its name
                                 //airLogic,                   //its mother
                                 //false,                         //no boulean operat
                                 //0);


 std::vector<G4VisAttributes*> fVisAttributes;
 auto visAttributes = new G4VisAttributes(G4Colour(0.,0.,0.)); // world is black
 //visAttributes->SetVisibility(false);
 fLogicWorld->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);                                     

 //visAttributes = new G4VisAttributes(G4Colour(0.4,0.4,0.4));   // make inner LAr detector dark grey
 //fLogicAbsorber->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 //visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,0.0));   // make Aluminium Layer yellow
 //fLogicAluminium->SetVisAttributes(visAttributes);
 //fVisAttributes.push_back(visAttributes);

 visAttributes = new G4VisAttributes(G4Colour(0.4,0.4,0.4,0.2));   // make outer LAr detector dark grey and opaque
 biggerLogicAbsorber->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

 visAttributes = new G4VisAttributes(G4Colour(0.19,0.83,0.78));   // make liquid scintilator turquoise
 lvol_fiber->SetVisAttributes(visAttributes);
 fVisAttributes.push_back(visAttributes);

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
    std::cout << pttoMaterial->GetName() << " has " << pttoMaterial->GetNumberOfMaterials() << " materials" << std::endl;
    G4MaterialTable * matTable = pttoMaterial->GetMaterialTable();
    G4int size = (int) (matTable->size());

    for (int i = 0 ; i < size ; i++){
        std::cout << (*matTable)[i] << std::endl;
    }

    const G4ElementVector * elmVect = pttoMaterial->GetElementVector();
    size = (int) (elmVect->size());

    for (int i = 0 ; i < size ; i++){
        std::cout << (*elmVect)[i] << std::endl;
    }

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

