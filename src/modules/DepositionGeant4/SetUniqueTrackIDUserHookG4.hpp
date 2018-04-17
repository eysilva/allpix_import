/**
 * @file
 * @brief Defines a user hook for Geant4 to assign custom track IDs which are unique
 * @copyright Copyright (c) 2018 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef SetUniqueTrackIDUserHookG4_H
#define SetUniqueTrackIDUserHookG4_H 1

#include "G4Track.hh"
#include "G4UserTrackingAction.hh"

#include "TrackInfoManager.hpp"

namespace allpix {
    /**
     * @brief Assigns every G4Track a TrackInfoG4 which carries the unique track id
     */
    class SetUniqueTrackIDUserHookG4 : public G4UserTrackingAction {
    public:
        /**
         * @brief Constructor taking a TrackInfoManager*
         * @param track_info_mgr_ptr Pointer to TrackInfoManager which must be used to create the TrackInfoG4 instances
         */
        SetUniqueTrackIDUserHookG4(TrackInfoManager* track_info_mgr_ptr) : track_info_mgr_ptr_(track_info_mgr_ptr){};

        /**
         * @brief Default destructor
         */
        ~SetUniqueTrackIDUserHookG4() = default;

        /**
         * @brief Default constructor is deleted
         */
        SetUniqueTrackIDUserHookG4() = delete;

        /**
         * @brief Called for every G4Track at beginning
         * @param aTrack The pointer to the G4Track for which this routine is called
         */
        void PreUserTrackingAction(const G4Track* aTrack);

        /**
         * @brief Called for every G4Track at end
         * @param aTrack The pointer to the G4Track for which this routine is called
         */
        void PostUserTrackingAction(const G4Track* aTrack);

    private:
        // Raw ptr to track info manager to create instances of TrackInfoG4
        TrackInfoManager* track_info_mgr_ptr_;
    };

} // namespace allpix
#endif
