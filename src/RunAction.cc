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
/// \file RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AccumulableManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    RunAction::RunAction()
    {
        // add new units for dose
        //
        const G4double milligray = 1.e-3 * gray;
        const G4double microgray = 1.e-6 * gray;
        const G4double nanogray = 1.e-9 * gray;
        const G4double picogray = 1.e-12 * gray;

        new G4UnitDefinition("milligray", "milliGy", "Dose", milligray);
        new G4UnitDefinition("microgray", "microGy", "Dose", microgray);
        new G4UnitDefinition("nanogray", "nanoGy", "Dose", nanogray);
        new G4UnitDefinition("picogray", "picoGy", "Dose", picogray);

        // Register accumulable to the accumulable manager
        G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
        accumulableManager->Register(fEdep);
        accumulableManager->Register(fEdep2);

        // Analysis manager
        auto* an = G4AnalysisManager::Instance();
        an->SetDefaultFileType("root");
        an->SetNtupleMerging(true);

        // H1 id 0: edep vs depth, 300 bins over 0..150nm (0.5nm spacing)
        an->CreateH1("EdepZ", "Energy deposition vs depth", 300, 0., 150.*nm);

        // Ntuple id 0: per-primary entry state at z = 0
        an->CreateNtuple("beam", "primary state at slab entry");
        an->CreateNtupleDColumn("x");    // col 0
        an->CreateNtupleDColumn("y");    // col 1
        an->CreateNtupleDColumn("xp");   // col 2, x' = px/pz
        an->CreateNtupleDColumn("yp");   // col 3
        an->CreateNtupleDColumn("E");    // col 4, kinetic energy
        an->FinishNtuple();
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::BeginOfRunAction(const G4Run*)
    {
        // inform the runManager to save random number seed
        G4RunManager::GetRunManager()->SetRandomNumberStore(false);

        // reset accumulables to their initial values
        G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
        accumulableManager->Reset();

        G4AnalysisManager::Instance()->OpenFile("target_run.root");

        // inside BeginOfRunAction(const G4Run*):
        if (IsMaster()) {
            G4EmCalculator calc;
            auto* ar40 = G4IonTable::GetIonTable()->GetIon(18, 40, 0.0);

            const std::vector<G4String> matNames = {"AmorphousCarbon", "G4_C"};
            const std::vector<G4double> energies = {
                25.*keV, 30.*keV, 35.*keV, 40.*keV, 41.*keV, 42.*keV, 43.*keV,
                44.*keV, 45.*keV, 46.*keV, 48.*keV, 50.*keV, 55.*keV, 60.*keV,
                70.*keV, 80.*keV, 100.*keV, 150.*keV, 200.*keV
            };

            for (const auto& mname : matNames) {
                auto* mat = G4Material::GetMaterial(mname);
                if (!mat) { G4cout << "Material " << mname << " not found\n"; continue; }
                G4cout << "\n=== dE/dx scan: Ar-40 in " << mname
                    << " (rho = " << mat->GetDensity()/(g/cm3) << " g/cm3) ===\n"
                    << "E[keV]   Se[keV/nm]   Sn[keV/nm]   total[keV/nm]\n";
                for (G4double e : energies) {
                    G4double se = calc.ComputeElectronicDEDX(e, ar40, mat);
                    G4double sn = calc.ComputeNuclearDEDX(e, ar40, mat);
                    G4cout << e/keV << "   " << se/(keV/nm) << "   "
                        << sn/(keV/nm) << "   " << (se+sn)/(keV/nm) << G4endl;
                }
            }
        }
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::EndOfRunAction(const G4Run* run)
    {
        G4int nofEvents = run->GetNumberOfEvent();
        if (nofEvents == 0) return;

        // Merge accumulables
        G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
        accumulableManager->Merge();

        // Write-out and close the analysis file
        auto* an = G4AnalysisManager::Instance();
        an->Write();
        an->CloseFile();

        // Compute dose = total energy deposit in a run and its variance
        //
        G4double edep = fEdep.GetValue();
        G4double edep2 = fEdep2.GetValue();

        G4double rms = edep2 - edep * edep / nofEvents;
        if (rms > 0.)
            rms = std::sqrt(rms);
        else
            rms = 0.;

        const auto detConstruction = static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        G4double mass = detConstruction->GetScoringVolume()->GetMass();
        G4double dose = edep / mass;
        G4double rmsDose = rms / mass;

        // Run conditions
        //  note: There is no primary generator action object for "master"
        //        run manager for multi-threaded mode.
        const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
            G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
        G4String runCondition;
        if (generatorAction) {
            const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
            runCondition += particleGun->GetParticleDefinition()->GetParticleName();
            runCondition += " of ";
            G4double particleEnergy = particleGun->GetParticleEnergy();
            runCondition += G4BestUnit(particleEnergy, "Energy");
        }

        // Print
        //
        if (IsMaster()) {
            G4cout << G4endl << "--------------------End of Global Run-----------------------";
        }
        else {
            G4cout << G4endl << "--------------------End of Local Run------------------------";
        }

        G4cout << G4endl << " The run is " << nofEvents << " " << runCondition << G4endl << G4endl;
        G4cout << "  --> cumulated edep per run in scoring volume = " << G4BestUnit(edep, "Energy") 
                << " = " << edep/joule << " joule" << G4endl;  
        G4cout << "  --> mass of scoring volume = " << G4BestUnit(mass, "Mass") << G4endl << G4endl; 
        G4cout << " Absorbed dose per run in scoring volume = edep/mass = " << G4BestUnit(dose, "Dose")
                << "; rms = " << G4BestUnit(rmsDose, "Dose") << G4endl
                << "------------------------------------------------------------" << G4endl << G4endl;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void RunAction::AddEdep(G4double edep)
    {
        fEdep += edep;
        fEdep2 += edep * edep;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
