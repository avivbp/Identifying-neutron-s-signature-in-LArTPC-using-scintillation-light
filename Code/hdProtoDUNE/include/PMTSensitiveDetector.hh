// PMTSensitiveDetector.hh
#include "G4VSensitiveDetector.hh"
class EventAction; // Forward declaration

class PMTSensitiveDetector : public G4VSensitiveDetector {
public:
    PMTSensitiveDetector(const G4String& name, EventAction* eventAction);
    ~PMTSensitiveDetector() override;

    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

private:
    EventAction* fEventAction; // Pointer to EventAction
};

