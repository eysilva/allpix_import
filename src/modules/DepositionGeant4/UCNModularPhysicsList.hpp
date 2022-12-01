#include "G4ProcessManager.hh"
#include "G4VPhysicsConstructor.hh"

namespace allpix {

    class UCNModularPhysicsList : public G4VPhysicsConstructor {

    public:
        UCNModularPhysicsList(const G4String& name = "ucn-physics");

        virtual ~UCNModularPhysicsList();
        virtual void ConstructParticle();
        virtual void ConstructProcess();
    };
}
