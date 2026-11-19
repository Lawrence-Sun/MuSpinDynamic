#pragma once

#include "../crmusr_common.hxx"
#include "TRandom3.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace stopped_muon {

// P = 2<S>/hbar in the laboratory frame; times in ns and fields in T.
enum class SpinState { MuPlus, MuMinus, Muonium };

struct SpinParameters {
    Vec3 external_field{0, 0, 0};
    double delta = 0;       // Local-field frequency width, ns^-1 (not cyclic MHz).
    double hop_rate = 0;    // ns^-1; zero selects the static model.
};

enum class DomainFieldMode { None, IsotropicDomains, SaturatedDomain };

// Static per-event internal magnetic field from ferromagnetic domains.
// Total static field seen by the muon is external_field + internal_field_T * m_hat,
// where internal_field_T is signed (negative when antiparallel to magnetization m_hat).
struct DomainFieldModel {
    DomainFieldMode mode = DomainFieldMode::None;
    double internal_field_T = 0;
    Vec3 magnetization_dir{0, 0, 1}; // Used only when mode == SaturatedDomain.
};

inline bool finite(const Vec3& v) {
    return std::isfinite(v.X()) && std::isfinite(v.Y()) && std::isfinite(v.Z());
}

inline void validate_domain_field(const DomainFieldModel& model) {
    if (!std::isfinite(model.internal_field_T)) {
        throw std::invalid_argument("domain internal_field_T must be finite");
    }
    if (model.mode == DomainFieldMode::SaturatedDomain) {
        if (!finite(model.magnetization_dir) || !std::isfinite(model.magnetization_dir.R()) ||
            model.magnetization_dir.R() <= 0) {
            throw std::invalid_argument("saturated domain magnetization_dir must be finite and non-zero");
        }
    } else if (model.mode != DomainFieldMode::None && model.mode != DomainFieldMode::IsotropicDomains) {
        throw std::invalid_argument("unknown domain field mode");
    }
}

inline Vec3 sample_domain_field(const DomainFieldModel& model, TRandom3& rng) {
    validate_domain_field(model);
    if (model.mode == DomainFieldMode::None || model.internal_field_T == 0) {
        return {0, 0, 0};
    }
    if (model.mode == DomainFieldMode::SaturatedDomain) {
        const Vec3 dir = model.magnetization_dir / model.magnetization_dir.R();
        const Vec3 field = model.internal_field_T * dir;
        if (!finite(field) || !std::isfinite(field.R())) {
            throw std::overflow_error("saturated domain field overflow");
        }
        return field;
    }
    // IsotropicDomains: draw a uniformly distributed unit vector m_hat on S^2.
    const double z = 2 * rng.Rndm() - 1;
    const double phi = 2 * pi * rng.Rndm();
    const double rho = std::sqrt(std::max(0., 1 - z * z));
    const Vec3 m_hat{rho * std::cos(phi), rho * std::sin(phi), z};
    const Vec3 field = model.internal_field_T * m_hat;
    if (!finite(field) || !std::isfinite(field.R())) {
        throw std::overflow_error("isotropic domain field overflow");
    }
    return field;
}

inline void require_nonnegative(double value, const char* name) {
    if (!std::isfinite(value) || value < 0) {
        throw std::invalid_argument(name);
    }
}

inline double frequency_scale(SpinState state) {
    switch (state) {
        case SpinState::MuPlus:
        case SpinState::MuMinus: return gamma_muon;
        case SpinState::Muonium: return gamma_muonium;
    }
    throw std::invalid_argument("unknown spin state");
}

inline double t1_relaxation_rate(double delta, double nu, double gamma, double b) {
    if (nu == 0) return 0;
    const double omega = gamma * b;
    return 2 * delta * delta * nu / (nu * nu + omega * omega);
}

inline double t2_relaxation_rate(double delta, double nu, double gamma, double b) {
    if (nu == 0) return 0;
    const double omega = gamma * b;
    return delta * delta / nu + delta * delta * nu / (nu * nu + omega * omega);
}

// These two helpers retain the legacy formula's +sin convention for regression.
inline Vec3 legacy_free_polarization(double t, double b, const Vec2& p) {
    const double phase = gamma_muon * b * t;
    return {p.X() * std::cos(phase), p.X() * std::sin(phase), p.Y()};
}

// Zeeman coefficients are coefficients of I_z in H/hbar (not magnetic moments).
inline Vec3 muonium_with_zeeman(double t, double b, const Vec2& p, double muon_zeeman) {
    const double h = 0.5 * omega_muonium_hfs;
    const double d = 0.5 * (gamma_electron - muon_zeeman) * b;
    const double s = 0.5 * (gamma_electron + muon_zeeman) * b;
    const double omega = std::sqrt(h * h + d * d);
    const double st = std::sin(s * t), ct = std::cos(s * t);
    const double so = std::sin(omega * t), co = std::cos(omega * t);
    const double transverse = p.X() * std::cos(h * t);
    return {transverse * (d * st * so + omega * ct * co) / omega,
            transverse * (omega * st * co - d * ct * so) / omega,
            p.Y() * (d * d + h * h * std::cos(2 * omega * t) + omega * omega)
                / (2 * omega * omega)};
}

inline Vec3 legacy_muonium_polarization(double t, double b, const Vec2& p) {
    return muonium_with_zeeman(t, b, p, gamma_muon);
}

// Prepared once per stopped particle. Evaluation has no RNG and is order-independent.
class SpinEvolution {
public:
    SpinEvolution(SpinState state, const Vec3& initial_polarization,
                  const SpinParameters& parameters, const Vec3& local_field = {})
        : state_(state), p0_(initial_polarization), local_field_(local_field) {
        const double gamma = frequency_scale(state);
        if (!finite(p0_) || !std::isfinite(p0_.R()) || p0_.R() > 1) {
            throw std::invalid_argument("initial laboratory polarization must have magnitude <= 1");
        }
        if (!finite(parameters.external_field) || !finite(local_field_)) {
            throw std::invalid_argument("magnetic fields must be finite");
        }
        require_nonnegative(parameters.delta, "delta must be finite and >= 0");
        require_nonnegative(parameters.hop_rate, "hop rate must be finite and >= 0");
        const Vec3 total = parameters.external_field + local_field_;
        b_ = total.R();
        const double external_b = parameters.external_field.R();
        if (!std::isfinite(b_) || !std::isfinite(external_b)) {
            throw std::invalid_argument("magnetic-field magnitude overflow");
        }
        z_ = b_ > 0 ? total / b_ : Vec3(0, 0, 1);
        parallel_ = p0_.Dot(z_);
        const Vec3 transverse = p0_ - parallel_ * z_;
        perpendicular_ = transverse.R();
        if (perpendicular_ > 0) {
            x_ = transverse / perpendicular_;
        } else if (z_.Rho() > 0) {
            x_ = {z_.Y() / z_.Rho(), -z_.X() / z_.Rho(), 0};
        } else {
            x_ = {1, 0, 0};
        }
        y_ = z_.Cross(x_);
        lambda1_ = t1_relaxation_rate(parameters.delta, parameters.hop_rate, gamma, external_b);
        lambda2_ = t2_relaxation_rate(parameters.delta, parameters.hop_rate, gamma, external_b);
        if (!std::isfinite(lambda1_) || !std::isfinite(lambda2_)) {
            throw std::invalid_argument("relaxation-rate overflow");
        }
    }

    Vec3 polarization(double time) const {
        require_nonnegative(time, "time must be finite and >= 0");
        if (time == 0) return p0_;
        const Vec2 p{perpendicular_ * std::exp(-lambda2_ * time),
                     parallel_ * std::exp(-lambda1_ * time)};
        Vec3 result = state_ == SpinState::Muonium
            ? muonium_with_zeeman(time, b_, p, -gamma_muon)
            : legacy_free_polarization(time, b_, p);
        // dP/dt = gamma_signed P x B: positive muons rotate toward -y for B=+z.
        if (state_ == SpinState::MuPlus) result.SetY(-result.Y());
        const Vec3 lab = result.X() * x_ + result.Y() * y_ + result.Z() * z_;
        if (!finite(lab)) throw std::overflow_error("non-finite spin evolution");
        return lab;
    }

    const Vec3& local_field() const { return local_field_; }
    double longitudinal_rate() const { return lambda1_; }
    double transverse_rate() const { return lambda2_; }

private:
    SpinState state_;
    Vec3 p0_, local_field_, x_, y_, z_;
    double b_, parallel_, perpendicular_, lambda1_, lambda2_;
};

inline SpinEvolution sample_spin_evolution(SpinState state, const Vec3& p0,
                                         const SpinParameters& parameters, TRandom3& rng) {
    // Validate before consuming random numbers, including when delta is zero.
    const SpinEvolution validated(state, p0, parameters);
    if (parameters.delta == 0) return validated;
    const double sigma = parameters.delta / frequency_scale(state);
    const Vec3 field{rng.Gaus(0, sigma), rng.Gaus(0, sigma), rng.Gaus(0, sigma)};
    return SpinEvolution(state, p0, parameters, field);
}

// Strong-collision trajectory for ordinary muons; no phenomenological damping.
class HoppingTrajectory {
public:
    struct Segment {
        double start;
        Vec3 initial_polarization;
        Vec3 local_field;
    };

    HoppingTrajectory(SpinState state, const Vec3& p0, const SpinParameters& parameters,
                      double duration, TRandom3& rng, std::size_t max_hops = 100000)
        : external_(parameters.external_field), duration_(duration),
          signed_gamma_(state == SpinState::MuMinus ? -gamma_muon : gamma_muon) {
        if (state != SpinState::MuPlus && state != SpinState::MuMinus) {
            throw std::invalid_argument("hopping supports ordinary muons only, not muonium");
        }
        // Reuse laboratory polarization / field validation without evaluating old rates.
        (void)SpinEvolution(state, p0, {external_, 0, 0});
        require_nonnegative(duration, "duration must be finite and >= 0");
        require_nonnegative(parameters.delta, "delta must be finite and >= 0");
        require_nonnegative(parameters.hop_rate, "hop rate must be finite and >= 0");
        if (parameters.hop_rate * duration > static_cast<double>(max_hops)) {
            throw std::length_error("expected hop count exceeds max_hops");
        }
        const double sigma = parameters.delta / gamma_muon;
        if (!std::isfinite(sigma)) throw std::invalid_argument("local-field width overflow");
        const auto sample_field = [&]() -> Vec3 {
            if (sigma == 0) return {};
            const Vec3 field{rng.Gaus(0, sigma), rng.Gaus(0, sigma), rng.Gaus(0, sigma)};
            if (!finite(field) || !std::isfinite((external_ + field).R())) {
                throw std::overflow_error("sampled magnetic-field overflow");
            }
            return field;
        };
        segments_.push_back({0, p0, sample_field()});
        if (parameters.hop_rate == 0 || duration == 0) return;
        while (true) {
            // Unit exponential avoids overflow in the mean 1/nu for tiny rates.
            const double wait = rng.Exp(1.) / parameters.hop_rate;
            const Segment& previous = segments_.back();
            if (wait > duration_ - previous.start) break;
            const double next = previous.start + wait;
            if (!(next > previous.start)) {
                throw std::overflow_error("hop time cannot advance at floating-point precision");
            }
            if (segments_.size() - 1 >= max_hops) {
                throw std::length_error("realized hop count exceeds max_hops; trajectory not returned");
            }
            const Vec3 p = propagate(previous, next - previous.start);
            segments_.push_back({next, p, sample_field()});
            if (next == duration_) break;
        }
    }

    Vec3 polarization(double time) const {
        require_nonnegative(time, "query time must be finite and >= 0");
        if (time > duration_) throw std::invalid_argument("query exceeds prepared duration");
        const auto after = std::upper_bound(segments_.begin(), segments_.end(), time,
            [](double t, const Segment& segment) { return t < segment.start; });
        const Segment& segment = *(after - 1);
        return propagate(segment, time - segment.start);
    }

    const std::vector<Segment>& segments() const { return segments_; }
    std::size_t hop_count() const { return segments_.size() - 1; }
    double duration() const { return duration_; }

private:
    Vec3 propagate(const Segment& segment, double elapsed) const {
        const Vec3 total = external_ + segment.local_field;
        const double b = total.R();
        if (elapsed == 0 || b == 0) return segment.initial_polarization;
        const Vec3 axis = total / b;
        const Vec3 parallel = segment.initial_polarization.Dot(axis) * axis;
        const Vec3 perpendicular = segment.initial_polarization - parallel;
        const double phase = signed_gamma_ * b * elapsed;
        const Vec3 result = parallel + perpendicular * std::cos(phase)
            + perpendicular.Cross(axis) * std::sin(phase);
        if (!finite(result)) throw std::overflow_error("non-finite hopping spin evolution");
        return result;
    }

    Vec3 external_;
    double duration_, signed_gamma_;
    std::vector<Segment> segments_;
};

} // namespace stopped_muon
