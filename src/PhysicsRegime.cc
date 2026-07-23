#include "PhysicsRegime.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4EmParameters.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4GenericIon.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4PhysicalConstants.hh"

namespace KACST 
{

    class IonLowEnergyFix  : public G4VPhysicsConstructor 
    {
    public:
        IonLowEnergyFix () : G4VPhysicsConstructor("IonLowEnergyFix") {}
        void ConstructParticle() override {}
        void ConstructProcess() override {
            auto* cfg = G4LossTableManager::Instance()->EmConfigurator();
            // 2 MeV/u boundary, mass-scaled to total kinetic energy (TestEm7 pattern)
            G4double eth = 2.*MeV * G4GenericIon::GenericIon()->GetPDGMass() / proton_mass_c2;
            cfg->SetExtraEmModel("GenericIon", "ionIoni",
                                new G4IonParametrisedLossModel(), "", 0.0, eth,
                                new G4IonFluctuations());
            cfg->SetExtraEmModel("GenericIon", "ionIoni",
                                new G4BetheBlochModel(), "", eth, 100*TeV,
                                new G4UniversalFluctuation());
        }
    };

    PhysicsRegime::PhysicsRegime()
    {
        // The one registered constructor: most accurate EM, includes
        // ICRU-table ion ionisation (Se) and G4NuclearStopping (Sn).
        RegisterPhysics(new G4EmStandardPhysics_option4());
        RegisterPhysics(new G4EmStandardPhysics_option3());

        RegisterPhysics(new G4StepLimiterPhysics()); // Step limiter to ensure step size is limited to the cut value, otherwise Geant4 will take large steps and miss energy deposition in the thin target.
        RegisterPhysics(new IonLowEnergyFix ()); 
        // No hadronic physics: at 25-50 keV there are no nuclear
        // reactions. Omission is correct, not an approximation.

        // Scale settings for a ~50 nm interaction depth.
        auto* em = G4EmParameters::Instance();
        // Energy-step per sim-iteration: 10% of range, 1 nm final range. Low-energy regime -> lower sim-steps or sim will be unstable
        em->SetStepFunctionIons(0.1, 1*nm);   // max 10% of range per step, 1 nm final range
        em->SetMinEnergy(100*eV);   // default Geant4 tracking floor is 1KeV, track ions below the 1 keV default floor

        // Depth profile-max ~100nm (carbon), gamme/energy deposition step-range default to 0.7mm(!), change to 1nm to ensure energy deposition is recorded in the 100nm target.
        SetDefaultCutValue(0.5*nm); // production range cut for secondaries
    }

    void PhysicsRegime::SetCuts()
    {
        // Allow the nm range cut to convert to energies below the 990 eV floor.
        G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 100*MeV);
        G4VModularPhysicsList::SetCuts();
    }
}