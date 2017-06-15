/**
 * @file
 * @brief Parameters of a pixel detector model
 *
 * @copyright MIT License
 */

#ifndef ALLPIX_PIXEL_DETECTOR_H
#define ALLPIX_PIXEL_DETECTOR_H

#include <string>
#include <utility>

#include <Math/Cartesian2D.h>
#include <Math/DisplacementVector2D.h>
#include <Math/Point3D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>

#include "DetectorModel.hpp"

// TODO [doc] This class is fully documented later when it is cleaned up further

namespace allpix {

    /**
     * @ingroup DetectorModels
     * @brief Model of a hybrid pixel detector: a detector where the sensor is bump-bonded to the chip
     */
    class HybridPixelDetectorModel : public DetectorModel {
    public:
        // Constructor and destructor
        explicit HybridPixelDetectorModel(std::string type)
            : DetectorModel(std::move(type)), number_of_pixels_(1, 1), coverlayer_material_("Al"), has_coverlayer_(false) {}

        /* Coordinate definitions
         * NOTE: center is at the middle of the first pixel at half z
         */
        // sensor coordinates relative to center of local frame
        double getSensorMinX() const override { return -getPixelSize().x() / 2.0; }
        double getSensorMinY() const override { return -getPixelSize().x() / 2.0; }
        double getSensorMinZ() const override { return -getSensorSize().z() / 2.0; }
        ROOT::Math::XYZPoint getCenter() const override {
            return ROOT::Math::XYZPoint(getSensorSize().x() / 2.0 - getPixelSize().x() / 2.0,
                                        getSensorSize().y() / 2.0 - getPixelSize().y() / 2.0,
                                        0);
        }

        /* Number of pixels */
        ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<int>> getNPixels() const override {
            return number_of_pixels_;
        }
        void setNPixels(ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<int>> val) {
            number_of_pixels_ = std::move(val);
        }

        /* Pixel dimensions */
        ROOT::Math::XYVector getPixelSize() const override { return pixel_size_; }
        void setPixelSize(ROOT::Math::XYVector val) { pixel_size_ = std::move(val); }

        /* Sensor offset */
        // FIXME: a PCB offset makes probably far more sense
        ROOT::Math::XYVector getSensorOffset() { return sensor_offset_; }
        double getSensorOffsetX() const { return sensor_offset_.x(); };
        double getSensorOffsetY() const { return sensor_offset_.y(); };
        double getSensorOffsetZ() const { return getHalfPCBSizeZ(); }; // FIXME: see relation with GetHalfWrapperDZ()

        void setSensorOffset(ROOT::Math::XYVector val) { sensor_offset_ = std::move(val); }

        /* Chip dimensions */
        ROOT::Math::XYZVector getChipSize() { return chip_size_; }
        double getChipSizeX() const { return chip_size_.x(); };
        double getChipSizeY() const { return chip_size_.y(); };
        double getChipSizeZ() const { return chip_size_.z(); };

        double getHalfChipSizeX() const { return chip_size_.x() / 2.0; };
        double getHalfChipSizeY() const { return chip_size_.y() / 2.0; };
        double getHalfChipSizeZ() const { return chip_size_.z() / 2.0; };

        void setChipSize(ROOT::Math::XYZVector val) { chip_size_ = std::move(val); }

        /* Chip offset */
        ROOT::Math::XYZVector getChipOffset() const { return chip_offset_; }
        double getChipOffsetX() const { return chip_offset_.x(); };
        double getChipOffsetY() const { return chip_offset_.y(); };
        double getChipOffsetZ() const { return chip_offset_.z(); };

        void setChipOffset(ROOT::Math::XYZVector val) { chip_offset_ = std::move(val); }

        /* PCB dimensions */
        ROOT::Math::XYZVector getPCBSize() { return pcb_size_; }
        double getPCBSizeX() const { return pcb_size_.x(); }
        double getPCBSizeY() const { return pcb_size_.y(); }
        double getPCBSizeZ() const { return pcb_size_.z(); }

        double getHalfPCBSizeX() const { return pcb_size_.x() / 2.0; };
        double getHalfPCBSizeY() const { return pcb_size_.y() / 2.0; };
        double getHalfPCBSizeZ() const { return pcb_size_.z() / 2.0; };

        void setPCBSize(ROOT::Math::XYZVector val) { pcb_size_ = std::move(val); }

        /* Bump bonds */
        double getBumpSphereRadius() const { return bump_sphere_radius_; }
        void setBumpSphereRadius(double val) { bump_sphere_radius_ = val; }

        double getBumpHeight() const { return bump_height_; }
        void setBumpHeight(double val) { bump_height_ = val; }

        ROOT::Math::XYVector getBumpOffset() const { return bump_offset_; }
        double getBumpOffsetX() const { return bump_offset_.x(); }
        double getBumpOffsetY() const { return bump_offset_.y(); }

        void setBumpOffset(ROOT::Math::XYVector val) { bump_offset_ = std::move(val); }

        double getBumpCylinderRadius() const { return bump_cylinder_radius_; }
        void setBumpCylinderRadius(double val) { bump_cylinder_radius_ = val; }

        /* Guard rings */
        double getGuardRingExcessTop() const { return guard_ring_excess_top_; }
        double getGuardRingExcessBottom() const { return guard_ring_excess_bottom_; }
        double getGuardRingExcessRight() const { return guard_ring_excess_right_; }
        double getGuardRingExcessLeft() const { return guard_ring_excess_left_; }

        void setGuardRingExcessTop(double val) { guard_ring_excess_top_ = val; }
        void setGuardRingExcessBottom(double val) { guard_ring_excess_bottom_ = val; }
        void setGuardRingExcessHRight(double val) { guard_ring_excess_right_ = val; }
        void getGuardRingExcessLeft(double val) { guard_ring_excess_left_ = val; }

        /* Wrapper calculations (FIXME: this has to be reworked...) */
        double getHalfWrapperDX() const { return getHalfPCBSizeX(); }
        double getHalfWrapperDY() const { return getHalfPCBSizeY(); }
        double getHalfWrapperDZ() const {

            double whdz = getHalfPCBSizeZ() + getHalfChipSizeZ() + getBumpHeight() / 2.0 + getSensorSize().z() / 2.0;

            if(has_coverlayer_) {
                whdz += getHalfCoverlayerHeight();
            }

            return whdz;
        }

        /* Coverlayer */
        bool hasCoverlayer() const { return has_coverlayer_; }
        double getCoverlayerHeight() const { return coverlayer_height_; }
        double getHalfCoverlayerHeight() const { return coverlayer_height_ / 2.0; }
        void setCoverlayerHeight(double val) {
            coverlayer_height_ = val;
            has_coverlayer_ = true;
        }

        std::string getCoverlayerMaterial() const { return coverlayer_material_; }
        void setCoverlayerMaterial(std::string mat) { coverlayer_material_ = std::move(mat); }

    private:
        ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<int>> number_of_pixels_;

        ROOT::Math::XYVector pixel_size_;

        ROOT::Math::XYVector sensor_offset_;

        ROOT::Math::XYZVector chip_size_;
        ROOT::Math::XYZVector chip_offset_;

        ROOT::Math::XYZVector pcb_size_;

        double bump_sphere_radius_{};
        double bump_height_{};
        ROOT::Math::XYVector bump_offset_;
        double bump_cylinder_radius_{};

        // FIXME: use a 4D vector here?
        double guard_ring_excess_top_{};
        double guard_ring_excess_bottom_{};
        double guard_ring_excess_right_{};
        double guard_ring_excess_left_{};

        double coverlayer_height_{};
        std::string coverlayer_material_;
        bool has_coverlayer_;
    };
} // namespace allpix

#endif /* ALLPIX_PIXEL_DETECTOR_H */
