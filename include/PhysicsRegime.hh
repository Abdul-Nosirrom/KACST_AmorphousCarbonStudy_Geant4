// **** PHYSICS REGIME DEFINITION DEFINING PROCESSES ****
//
// Low-Energy Regime, Nuclear-Stopping Dominant (m=1/2)
// - Electronic Stopping: Via G4ionIonisation
// - Nuclear Stopping: Via G4NuclearStopping
//
// Default QBBC invalid, high-energy hadronic physics. Ensure no active hadronic processes are registered in this regime.



#pragma once
#include "G4VModularPhysicsList.hh"

namespace KACST
{
    class PhysicsRegime : public G4VModularPhysicsList
    {
        public:
            PhysicsRegime();
            ~PhysicsRegime() override = default;

            void SetCuts() override;
    };
}