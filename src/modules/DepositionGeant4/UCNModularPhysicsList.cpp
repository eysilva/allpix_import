/*
 * UCNPhysics.cpp
 *
 *  Created on: May 20, 2019
 *      Author: Mathieu Benoit
 */

#include "UCNModularPhysicsList.hpp"
#include "G4UCNAbsorption.hh"
#include "G4UCNBoundaryProcess.hh"
#include "G4UCNLoss.hh"
#include "G4UCNMultiScattering.hh"

using namespace allpix;

UCNModularPhysicsList::UCNModularPhysicsList(const G4String& name) : G4VPhysicsConstructor(name) {}

UCNModularPhysicsList::~UCNModularPhysicsList() {}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

void UCNModularPhysicsList::ConstructParticle() {}

void UCNModularPhysicsList::ConstructProcess() {

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    G4ProcessManager* pmanager = nullptr;

    while((*particleIterator)()) {
        G4ParticleDefinition* particle = particleIterator->value();
        pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if(!pmanager) {
            std::ostringstream o;
            o << "Particle " << particleName << "without a Process Manager";
            G4Exception("ExUCNExtraPhysics::ConstructProcess()", "", FatalException, o.str().c_str());
        }

        if(particleName == "neutron") {
            pmanager->AddDiscreteProcess(new G4UCNLoss());
            pmanager->AddDiscreteProcess(new G4UCNAbsorption());
            pmanager->AddDiscreteProcess(new G4UCNMultiScattering());

            G4UCNBoundaryProcess* ucnBoundaryProcess = new G4UCNBoundaryProcess();
            ucnBoundaryProcess->SetMicroRoughness(true);
            ucnBoundaryProcess->SetVerboseLevel(0);

            pmanager->AddDiscreteProcess(ucnBoundaryProcess);
        }
    }
}
