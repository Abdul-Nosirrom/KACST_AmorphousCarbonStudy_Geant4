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
/// \file SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void SteppingAction::UserSteppingAction(const G4Step* step)
    {
        if (!fScoringVolume) {
            const auto det = static_cast<const DetectorConstruction*>(
                G4RunManager::GetRunManager()->GetUserDetectorConstruction());
            fScoringVolume = det->GetScoringVolume();
        }

        auto* an = G4AnalysisManager::Instance();

        // record primary's endpoint: fires on the step where the primary dies
        auto* trk = step->GetTrack();
        if (trk->GetParentID() == 0 &&
            trk->GetTrackStatus() != fAlive) {
            auto p = step->GetPostStepPoint()->GetPosition();
            fEventAction->SetPrimaryEndPos(p.x(), p.y(), p.z());
}
        // (2) FIRST: beam state at slab entry -- must run before any
        // scoring-volume early return, because the entry step's PRE volume
        // is the World, not the slab.
        auto* post = step->GetPostStepPoint();
        if (step->GetTrack()->GetParentID() == 0 &&
            post->GetStepStatus() == fGeomBoundary &&
            post->GetTouchableHandle()->GetVolume() &&
            post->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fScoringVolume &&
            step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() != fScoringVolume) {
            auto pos = post->GetPosition();
            auto dir = post->GetMomentumDirection();
            an->FillNtupleDColumn(0, 0, pos.x()/mm);
            an->FillNtupleDColumn(0, 1, pos.y()/mm);
            an->FillNtupleDColumn(0, 2, dir.x()/dir.z());
            an->FillNtupleDColumn(0, 3, dir.y()/dir.z());
            an->FillNtupleDColumn(0, 4, post->GetKineticEnergy()/keV);
            an->AddNtupleRow(0);
        }

        // NOW the scoring-volume gate for everything edep-related
        G4LogicalVolume* volume =
            step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        if (volume != fScoringVolume) return;

        G4double edep = step->GetTotalEnergyDeposit();
        G4double niel = step->GetNonIonizingEnergyDeposit();   // nuclear-channel share of this step
        fEventAction->AddEdep(edep);
        if (edep > 0.) {
            G4double z = 0.5*(step->GetPreStepPoint()->GetPosition().z()
                            + step->GetPostStepPoint()->GetPosition().z());
            an->FillH1(0, z, edep);
            if (niel > 0.)        an->FillH1(1, z, niel);
            if (edep - niel > 0.) an->FillH1(2, z, edep - niel);
        }
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
