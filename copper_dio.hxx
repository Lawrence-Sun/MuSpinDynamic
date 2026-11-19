#pragma once

#include "bound_dio.hxx"

#include <stdexcept>
#include <utility>

namespace stopped_muon {

namespace cu_dio_detail {
inline constexpr double electron_mass_MeV = bound_dio_detail::electron_mass_MeV;
} // namespace cu_dio_detail

// Backward-compatible Copper DIO spectrum wrapper delegating to BoundDioSpectrum.
// Certified shape only: Q never changes caller-supplied competing rates.
class CopperDioSpectrum {
public:
    explicit CopperDioSpectrum(int isotope)
        : inner_(resolve_copper_nuclide(isotope)) {}

    int isotope() const { return inner_.isotope(); }
    double endpoint_MeV() const { return inner_.endpoint_MeV(); }
    double Q_shape() const { return inner_.Q_shape(); }
    const char* table_sha256() const { return inner_.table_sha256(); }
    double envelope_integral() const { return inner_.envelope_integral(); }

    std::pair<double, double> density(double energy) const {
        return inner_.density(energy);
    }

    DioElectron sample(
            const Vec3& polarization, TRandom3& rng, std::size_t max_trials = 100000) const {
        return inner_.sample(polarization, rng, max_trials);
    }

private:
    static DioNuclide resolve_copper_nuclide(int isotope) {
        if (isotope == 63) return DioNuclide::Cu63;
        if (isotope == 65) return DioNuclide::Cu65;
        throw std::invalid_argument("Cu DIO isotope must be explicitly 63 or 65");
    }

    BoundDioSpectrum inner_;
};

} // namespace stopped_muon
