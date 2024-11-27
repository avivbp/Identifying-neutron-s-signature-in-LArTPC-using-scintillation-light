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
/// \file electromagnetic/TestEm5/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Colour.hh"
#include "G4Scintillation.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include <sstream>

class G4Box;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class G4UniformMagField;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  explicit DetectorConstruction();
  virtual ~DetectorConstruction();

  void SetAbsorberMaterial (const G4String&);
  void SetAbsorberThickness(G4double);
  void SetAbsorberSizeYZ   (G4double);

  void SetAbsorberXpos(G4double);

  void SetWorldMaterial(const G4String&);
  void SetWorldSizeX   (G4double);
  void SetWorldSizeYZ  (G4double);

  void SetMagField(G4double);

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  void PrintGeomParameters();

  const G4Material* GetAbsorberMaterial() const {return fAbsorberMaterial;};
  G4double GetAbsorberThickness() const         {return fAbsorberThickness;};
  G4double GetAbsorberSizeYZ() const            {return fAbsorberSizeYZ;};

  G4double GetAbsorberXpos() const              {return fXposAbs;};
  G4double GetxstartAbs() const                 {return fXstartAbs;};
  G4double GetxendAbs() const                   {return fXendAbs;};
  
  std::vector<G4double> absorptionLength;
  G4double birks;

  const G4Material* GetWorldMaterial() const    {return fWorldMaterial;};
  G4double GetWorldSizeX() const                {return fWorldSizeX;};

  const G4VPhysicalVolume* GetAbsorber() const  {return fPhysiAbsorber;};

  const G4ThreeVector & MatVectMul(G4RotationMatrix mat,const G4ThreeVector & vect){

      G4double x = mat.xx()*vect.x() + mat.xy()*vect.y()+mat.xz()*vect.z();
      G4double y = mat.yx()*vect.x() + mat.yy()*vect.y()+mat.yz()*vect.z();
      G4double z = mat.zx()*vect.x() + mat.zy()*vect.y()+mat.zz()*vect.z();

      const G4ThreeVector & ret = G4ThreeVector(x,y,z);
      std::cout << "ret = " << ret << std::endl;
      return ret;
  }

  const G4ThreeVector & rotateAroundY(const G4ThreeVector & v,G4double angle){

      G4ThreeVector x = G4ThreeVector(std::cos(angle), 0, std::sin(-angle));
      G4ThreeVector y = G4ThreeVector(0, 1, 0);
      G4ThreeVector z = G4ThreeVector(std::sin(angle),0,std::cos(angle));
      G4RotationMatrix rotm1  = G4RotationMatrix(x, y, z);

      const G4ThreeVector & ret = MatVectMul(rotm1,v);
      return ret;
  }
 
  const void placeLiquidScintillator(G4ThreeVector pos,G4double rotatAng,G4double fiberDiam,G4double fiberLen,G4Material* fibMat,G4LogicalVolume* motherVol,G4String name){
      
      G4Tubs* solidFiber = new G4Tubs(name,                      //name
                             0*CLHEP::mm, 0.5*fiberDiam,       //r1, r2
                             0.5*fiberLen,               //half-length
                             0., CLHEP::twopi);                    //theta1, theta2

      G4LogicalVolume* logicFiber = new G4LogicalVolume(solidFiber,          //solid
                                       fibMat,            //material
                                       name);            //name


      // rotation around y axis
      // u, v, w are the daughter axes, projected on the mother frame
      G4double rotatAng2 = CLHEP::pi / 2 - rotatAng;
      G4ThreeVector u = G4ThreeVector(std::cos(rotatAng2) , 0, std::sin(rotatAng2));
      G4ThreeVector v = G4ThreeVector(0 , 1, 0);
      G4ThreeVector w = G4ThreeVector(std::sin(-rotatAng2), 0, std::cos(rotatAng2)); 


      // rotation around x axis
      //G4ThreeVector u1 = G4ThreeVector(1,0,0);
      //G4ThreeVector v1 = G4ThreeVector(0,std::cos(rotatAng),std::sin(-rotatAng));
      //G4ThreeVector w1 = G4ThreeVector(0,std::sin(rotatAng),std::cos(rotatAng));

      G4RotationMatrix rotMat  = G4RotationMatrix(u, v, w);
      std::cout << "position before rotation = " << pos << std::endl;
      pos = rotateAroundY(pos,rotatAng);
      std::cout << "position after rotation = " << pos << std::endl;
      const G4Transform3D & trans = G4Transform3D(rotMat,pos);


      G4VPhysicalVolume* physiFiber = new G4PVPlacement(trans,
                                    logicFiber,                    //its logical volume
                                    name,                          //its name
                                    motherVol,                      //its mother
                                    false,                         //no boulean operat
                                    0);      
  }

private:

  void DefineMaterials();
  void ComputeGeomParameters();
  void ChangeGeometry();

  G4double           fiberDiameter;
  G4double           fiberLength;
  G4Material *       fiberMat;
  G4Material*        Al;

  G4Material*        fAbsorberMaterial;
  G4double           fAbsorberThickness;
  G4double           fAbsorberSizeYZ;

  G4double           fXposAbs;
  G4double           fXstartAbs, fXendAbs;

  G4Material*        fWorldMaterial;
  G4double           fWorldSizeX;
  G4double           fWorldSizeYZ;

  G4Box*             fSolidWorld;
  G4LogicalVolume*   fLogicWorld;
  G4VPhysicalVolume* fPhysiWorld;

  G4Box*             airSolid;
  G4LogicalVolume*   airLogic;
  G4VPhysicalVolume* airPhysi;

  G4Box*             hydroSolid;
  G4LogicalVolume*   hydroLogic;
  G4VPhysicalVolume* hydroPhysi;

  G4Orb*            biggerSolidAbsorber;
  G4LogicalVolume*   biggerLogicAbsorber;
  G4VPhysicalVolume* biggerPhysiAbsorber;

  G4Tubs*            fSolidTeflon;
  G4LogicalVolume*   fLogicTeflon;
  G4VPhysicalVolume* fPhysiTeflon;

 // G4Box*             fSolidAbsorber;
  G4Tubs*            fSolidAbsorber;
  G4LogicalVolume*   fLogicAbsorber;
  G4VPhysicalVolume* fPhysiAbsorber;
 
  G4Tubs*          svol_fiber;
  G4LogicalVolume* lvol_fiber;
  G4VPhysicalVolume* pvol_fiber;
 
  G4Tubs*          new_svol_fiber;
  G4LogicalVolume* new_lvol_fiber;
  G4VPhysicalVolume* new_pvol_fiber;

   
  DetectorMessenger* fDetectorMessenger;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

