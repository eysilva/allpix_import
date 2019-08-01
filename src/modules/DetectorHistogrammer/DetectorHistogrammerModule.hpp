/**
 * @file
 * @brief Definition of detector histogramming module
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_MODULE_DETECTOR_HISTOGRAMMER_H
#define ALLPIX_MODULE_DETECTOR_HISTOGRAMMER_H

#include <memory>

#include <string>
#include <vector>

#include <ROOT/TThreadedObject.hxx>
#include <TH1I.h>
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include "core/config/Configuration.hpp"
#include "core/geometry/GeometryManager.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/module/Event.hpp"
#include "core/module/Module.hpp"

#include "Cluster.hpp"
#include "objects/PixelHit.hpp"

namespace allpix {
    /**
     * @ingroup Modules
     * @brief Module to plot the final digitized pixel data
     *
     * Generates a hitmap of all the produced pixel hits, together with a histogram of the cluster size
     */
    class DetectorHistogrammerModule : public Module {
    public:
        /**
         * @brief Constructor for this detector-specific module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param messenger Pointer to the messenger object to allow binding to messages on the bus
         * @param detector Pointer to the detector for this module instance
         */
        DetectorHistogrammerModule(Configuration& config, Messenger*, std::shared_ptr<Detector>);

        /**
         * @brief Initialize the histograms
         */
        void init(std::mt19937_64&) override;

        /**
         * @brief Fill the histograms
         */
        void run(Event*) override;

        /**
         * @brief Write the histograms to the modules file
         */
        void finalize() override;

    private:
        /**
         * @brief Perform a sparse clustering on the PixelHits
         */
        static std::vector<Cluster> doClustering(std::shared_ptr<PixelHitMessage>& pixels_message);

        /**
         * @brief analyze the available MCParticles and return the all particles identified as primary (i.e. that do not have
         * a parent). This might be several particles.
         */
        static std::vector<const MCParticle*> getPrimaryParticles(std::shared_ptr<MCParticleMessage>& mcparticle_message);

        std::shared_ptr<Detector> detector_;

        // Statistics to compute mean position
        ROOT::Math::XYVector total_vector_{};
        unsigned long total_hits_{};

        // Cut criteria for efficiency measurement:
        ROOT::Math::XYVector matching_cut_{};

        // Reference track resolution
        ROOT::Math::XYVector track_resolution_{};

        // Histograms to output
        ROOT::TThreadedObject<TH2D>*hit_map, *charge_map, *cluster_map;
        ROOT::TThreadedObject<TProfile2D>*cluster_size_map, *cluster_size_x_map, *cluster_size_y_map;
        ROOT::TThreadedObject<TProfile2D>*cluster_charge_map, *seed_charge_map;
        ROOT::TThreadedObject<TProfile2D>*residual_map, *residual_x_map, *residual_y_map;
        ROOT::TThreadedObject<TH1D>*residual_x, *residual_y;
        ROOT::TThreadedObject<TProfile>*residual_x_vs_x, *residual_y_vs_y, *residual_x_vs_y, *residual_y_vs_x;
        ROOT::TThreadedObject<TProfile2D>*efficiency_map, *efficiency_detector;
        ROOT::TThreadedObject<TProfile>*efficiency_vs_x, *efficiency_vs_y;
        ROOT::TThreadedObject<TH1D>* event_size;
        ROOT::TThreadedObject<TH1D>* cluster_size;
        ROOT::TThreadedObject<TH1D>* cluster_size_x;
        ROOT::TThreadedObject<TH1D>* cluster_size_y;
        ROOT::TThreadedObject<TH1D>* n_cluster;
        ROOT::TThreadedObject<TH1D>* cluster_charge;
    };
} // namespace allpix

#endif /* ALLPIX_MODULE_DETECTOR_HISTOGRAMMER_H */
