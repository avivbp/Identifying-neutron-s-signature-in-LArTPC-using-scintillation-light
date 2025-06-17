#include "PMTSensitiveDetector.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"

PMTSensitiveDetector::PMTSensitiveDetector(const G4String& name, EventAction* eventAction)
    : G4VSensitiveDetector(name), fEventAction(eventAction) {}

PMTSensitiveDetector::~PMTSensitiveDetector() = default;

G4bool PMTSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {

    G4cout << "ProcessHits called for step in sensitive detector." << G4endl;
    // Get the track associated with this step
    G4Track* track = step->GetTrack();

    // Check if the particle is an optical photon
    if (track->GetDefinition() == G4OpticalPhoton::Definition()) {
        // Increment the photon count in EventAction
        if (fEventAction) {
            std::cout << "hit PMT" << std::endl;
            fEventAction->num += 1;
        } else {
            G4cerr << "Error: EventAction pointer is null in PMTSensitiveDetector!" << G4endl;
        }

        // Optionally, terminate the photon to avoid unnecessary tracking
        // track->SetTrackStatus(fStopAndKill);

        return true; // Signal that the hit was processed successfully
    }

    return false; // Not an optical photon; no processing done
}

