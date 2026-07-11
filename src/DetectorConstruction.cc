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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include <G4Types.hh>

namespace B1 
{
    G4VPhysicalVolume* DetectorConstruction::Construct() 
    {
        // Get nist material manager
        G4NistManager *nist = G4NistManager::Instance();

        // Option to switch on/off checking of volumes overlaps. In-case.
        G4bool checkOverlaps = true;

        // World: Vacuum in cm, chamber size or bigger than beam + target dimensions
        G4double world_sizeXY = 20 * cm, world_sizeZ = 30 * cm; // Size arbitrary, but too small is bad so need to be big enough
        G4Material *vacuum_mat = nist->FindOrBuildMaterial("G4_Galactic");
        G4Box* world_bounds = new G4Box("World", 0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ); // Box half-extents

        auto world_log = new G4LogicalVolume(world_bounds, vacuum_mat, "World");
        // Place it
        auto physWorld = new G4PVPlacement(0, G4ThreeVector(), world_log, "World", 0, false, 0, checkOverlaps);

        // Carbon-Slab, amorphous carbon w/ denisty 1.8-2.2 g/cm3, 1.9 g/cm3 is used here based on SRIM test + Gupta paper
        // TODO: Setup messanger to make it a variable for density, to sweep the 1.8-2.2 range
        G4Material* carbon_mat = nist->BuildMaterialWithNewDensity("AmorphousCarbon", "G4_C", 1.9 * g/cm3); // 1.9 g/cm3 is used here based on SRIM test + Gupta paper
        G4double carbon_thickness = 150 * nm; // From SRIM range profile, roughly 34nm mean, and low straggling, so 150nm is acceptable and matches roughly gupta paper
        G4double carbon_size_XY = 10 * cm; // roughly
        G4double carbon_density = 1.9 * g / cm3;
        G4Box* carbon_slab_dim = new G4Box("CarbonSlab", 0.5 * carbon_size_XY, 0.5 * carbon_size_XY, 0.5 * carbon_thickness);
        // Placement - positions box-center, Z-face of carbon slab should be at Z=0, so center is at 0.5*thickness
        G4ThreeVector carbon_slab_pos = G4ThreeVector(0, 0, 0.5 * carbon_thickness); // Beam travels +Z, placement has Z = 0 be the front-face
        auto carbon_slab_log = new G4LogicalVolume(carbon_slab_dim, carbon_mat, "CarbonSlab");
        new G4PVPlacement(0, carbon_slab_pos, carbon_slab_log, "CarbonSlab", world_log, false, 0, checkOverlaps);

        carbon_slab_log->SetUserLimits(new G4UserLimits( 0.5 * nm)); // 0.5nm steps through the carbon volume

        // Carbon slab is the scoring volume
        fScoringVolume = carbon_slab_log;

        return physWorld;
    }

} // namespace B1
