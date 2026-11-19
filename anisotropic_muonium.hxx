#pragma once

#include "muonium_hopping.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace stopped_muon {

namespace anisotropic_muonium_detail {

inline Vec3 sample_unit_sphere(TRandom3& rng) {
    const double z = 2 * rng.Rndm() - 1;
    const double phi = 2 * pi * rng.Rndm();
    const double rho = std::sqrt(std::max(0., 1 - z * z));
    return {rho * std::cos(phi), rho * std::sin(phi), z};
}

// Trace-preserving, completely positive electron-spin randomization:
//   rho -> Tr_e(rho) (x) (I_e / 2)
// Preserves Tr(rho) = 1, Hermiticity, positive semidefiniteness, and exact muon polarization P_mu.
inline DensityMatrix randomize_electron_spin(const DensityMatrix& rho) {
    DensityMatrix out{};
    for (int a = 0; a < 2; ++a) {
        for (int ap = 0; ap < 2; ++ap) {
            const std::complex<double> rho_mu =
                rho[2 * a][2 * ap] + rho[2 * a + 1][2 * ap + 1];
            out[2 * a][2 * ap] = 0.5 * rho_mu;
            out[2 * a + 1][2 * ap + 1] = 0.5 * rho_mu;
        }
    }
    return out;
}

// Axial hyperfine + Zeeman Hamiltonian in rad/ns:
//   H / hbar = A_iso (I_mu . S_e)
//            + D [ (I_mu . n_hat)(S_e . n_hat) - (1/3)(I_mu . S_e) ]
//            + B . (-gamma_muon I_mu + gamma_electron S_e)
// Principal hyperfine values:
//   A_parallel = A_iso + (2/3) D,   A_perp = A_iso - (1/3) D,   A_parallel - A_perp = D.
// Both A_iso and D are signed finite real values (rad/ns).
inline DensityMatrix axial_hamiltonian(
        const Vec3& field, double A_iso, double D, const Vec3& axis_unit) {
    if (!finite(field) || !std::isfinite(field.R())) {
        throw std::invalid_argument("anisotropic muonium magnetic field must be finite");
    }
    if (!std::isfinite(A_iso) || !std::isfinite(D)) {
        throw std::invalid_argument("anisotropic muonium A_iso and D must be finite");
    }
    if (!finite(axis_unit) || !std::isfinite(axis_unit.R()) || axis_unit.R() <= 0) {
        throw std::invalid_argument("anisotropic muonium symmetry axis must be finite and non-zero");
    }
    const Vec3 n = axis_unit / axis_unit.R();
    const double nv[3] = {n.X(), n.Y(), n.Z()};

    // Symmetric 3x3 hyperfine tensor T_ab = (A_iso - D/3) delta_ab + D n_a n_b
    double T[3][3]{};
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            T[a][b] = ((a == b) ? (A_iso - D / 3.0) : 0.0) + D * nv[a] * nv[b];
        }
    }

    // Spin-1/2 matrices S_x, S_y, S_z in {|up>, |down>} basis
    using C = std::complex<double>;
    const C S[3][2][2] = {
        {{C(0, 0), C(0.5, 0)}, {C(0.5, 0), C(0, 0)}},
        {{C(0, 0), C(0, -0.5)}, {C(0, 0.5), C(0, 0)}},
        {{C(0.5, 0), C(0, 0)}, {C(0, 0), C(-0.5, 0)}},
    };

    DensityMatrix h{};
    for (int am = 0; am < 2; ++am) {
        for (int ae = 0; ae < 2; ++ae) {
            const int i = 2 * am + ae;
            for (int bm = 0; bm < 2; ++bm) {
                for (int be = 0; be < 2; ++be) {
                    const int j = 2 * bm + be;
                    C val(0, 0);
                    for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            if (T[a][b] != 0.0) {
                                val += T[a][b] * S[a][am][bm] * S[b][ae][be];
                            }
                        }
                    }
                    h[i][j] = val;
                }
            }
        }
    }

    // Add Zeeman terms matching muonium_detail::hamiltonian exactly
    for (int i = 0; i < 4; ++i) {
        const double m = i < 2 ? 0.5 : -0.5;
        const double e = i % 2 == 0 ? 0.5 : -0.5;
        h[i][i] += field.Z() * (-gamma_muon * m + gamma_electron * e);
        h[i][i ^ 2] += -gamma_muon / 2.0 * C(field.X(), -2.0 * m * field.Y());
        h[i][i ^ 1] += gamma_electron / 2.0 * C(field.X(), -2.0 * e * field.Y());
    }
    return h;
}

} // namespace anisotropic_muonium_detail

struct AnisotropicMuoniumParameters {
    Vec3 external_field{0, 0, 0};     // T
    double local_field_sigma_T = 0;   // Per Cartesian component, T
    double hop_rate = 0;              // Local-field jump rate, ns^-1
    double hyperfine_iso = omega_muonium_hfs; // Signed isotropic A_iso in rad/ns
    double dipolar_D = 0;             // Signed axial dipolar D = A_parallel - A_perp in rad/ns
    double spin_exchange_rate = 0;    // Electron spin randomization rate, ns^-1
    double reorientation_rate = 0;    // Molecular symmetry axis jump rate, ns^-1
    Vec3 fixed_axis{0, 0, 0};         // Zero => random isotropic axis per event (powder average)
};

class AnisotropicMuoniumTrajectory;

// Single-segment coherent evolution under an axial hyperfine tensor + magnetic field.
class AnisotropicMuoniumEvolution {
public:
    AnisotropicMuoniumEvolution(
            const Vec3& p0, const Vec3& field, double A_iso, double D, const Vec3& axis)
        : AnisotropicMuoniumEvolution(muonium_detail::initial_density(p0), field, A_iso, D, axis) {}

    DensityMatrix density(double time) const {
        require_nonnegative(time, "anisotropic muonium time must be finite and >= 0");
        if (time == 0) return initial_;
        DensityMatrix u{};
        for (int k = 0; k < 8; ++k) {
            const double phase = energies_[k] * time;
            if (!std::isfinite(phase)) {
                throw std::overflow_error("anisotropic muonium phase overflow");
            }
            const double c = std::cos(phase), s = std::sin(phase);
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    const double a = eigenvectors_[i][k] * eigenvectors_[j][k];
                    const double b = eigenvectors_[i + 4][k] * eigenvectors_[j][k];
                    u[i][j] += std::complex<double>(a * c + b * s, b * c - a * s);
                }
            }
        }
        DensityMatrix temporary{}, result{};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    temporary[i][j] += u[i][k] * initial_[k][j];
                }
            }
        }
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    result[i][j] += temporary[i][k] * std::conj(u[j][k]);
                }
            }
        }
        (void)muonium_detail::real_embedding(result);
        return result;
    }

    Vec3 polarization(double time) const { return muon_polarization(density(time)); }
    const DensityMatrix& initial_density() const { return initial_; }

private:
    friend class AnisotropicMuoniumTrajectory;
    AnisotropicMuoniumEvolution(
            const DensityMatrix& initial, const Vec3& field, double A_iso, double D, const Vec3& axis)
        : initial_(initial) {
        const TMatrixDSymEigen spectrum(muonium_detail::real_embedding(
            anisotropic_muonium_detail::axial_hamiltonian(field, A_iso, D, axis)));
        const auto& values = spectrum.GetEigenValues();
        const auto& vectors = spectrum.GetEigenVectors();
        for (int k = 0; k < 8; ++k) {
            energies_[k] = values[k];
            if (!std::isfinite(energies_[k])) {
                throw std::overflow_error("anisotropic muonium spectrum overflow");
            }
            for (int i = 0; i < 8; ++i) {
                eigenvectors_[i][k] = vectors(i, k);
                if (!std::isfinite(eigenvectors_[i][k])) {
                    throw std::overflow_error("anisotropic muonium eigenvector overflow");
                }
            }
        }
    }

    DensityMatrix initial_;
    std::array<double, 8> energies_{};
    std::array<std::array<double, 8>, 8> eigenvectors_{};
};

// Stochastic 4x4 trajectory with axial hyperfine tensor, optional powder averaging,
// local-field hopping, strong-collision electron spin exchange, and molecular reorientation.
class AnisotropicMuoniumTrajectory {
public:
    struct Segment {
        double start;
        Vec3 local_field;
        Vec3 symmetry_axis;
        AnisotropicMuoniumEvolution evolution;
    };

    AnisotropicMuoniumTrajectory(
            const Vec3& p0, const AnisotropicMuoniumParameters& parameters,
            double duration, TRandom3& rng, std::size_t max_hops = 100000)
        : duration_(duration) {
        // Validate all deterministic inputs before consuming any RNG draws.
        require_nonnegative(duration, "anisotropic muonium duration must be finite and >= 0");
        require_nonnegative(parameters.local_field_sigma_T, "local field sigma must be finite and >= 0");
        require_nonnegative(parameters.hop_rate, "hop rate must be finite and >= 0");
        require_nonnegative(parameters.spin_exchange_rate, "spin exchange rate must be finite and >= 0");
        require_nonnegative(parameters.reorientation_rate, "reorientation rate must be finite and >= 0");
        if (!finite(parameters.fixed_axis) || !std::isfinite(parameters.fixed_axis.R())) {
            throw std::invalid_argument("anisotropic muonium fixed_axis must be finite");
        }
        const bool use_fixed_axis = parameters.fixed_axis.R() > 0;
        const Vec3 probe_axis = use_fixed_axis
            ? (parameters.fixed_axis / parameters.fixed_axis.R())
            : Vec3{0, 0, 1};
        const AnisotropicMuoniumEvolution validated(
            p0, parameters.external_field, parameters.hyperfine_iso, parameters.dipolar_D, probe_axis);
        const double total_rate =
            parameters.hop_rate + parameters.spin_exchange_rate + parameters.reorientation_rate;
        if (!std::isfinite(total_rate)) {
            throw std::invalid_argument("anisotropic muonium total jump rate overflow");
        }
        if (total_rate * duration > static_cast<double>(max_hops)) {
            throw std::length_error("expected anisotropic muonium jump count exceeds max_hops");
        }

        const auto sample_field = [&]() -> Vec3 {
            if (parameters.local_field_sigma_T == 0) return {0, 0, 0};
            const double s = parameters.local_field_sigma_T;
            const Vec3 f{rng.Gaus(0, s), rng.Gaus(0, s), rng.Gaus(0, s)};
            if (!finite(f) || !std::isfinite((parameters.external_field + f).R())) {
                throw std::overflow_error("sampled anisotropic muonium field overflow");
            }
            return f;
        };

        Vec3 current_axis = use_fixed_axis
            ? probe_axis
            : anisotropic_muonium_detail::sample_unit_sphere(rng);
        Vec3 current_field = sample_field();
        DensityMatrix rho = validated.initial_density();
        segments_.push_back({
            0.0, current_field, current_axis,
            AnisotropicMuoniumEvolution(
                rho, parameters.external_field + current_field,
                parameters.hyperfine_iso, parameters.dipolar_D, current_axis)
        });

        if (total_rate == 0 || duration == 0) return;

        while (true) {
            const double wait = rng.Exp(1.0) / total_rate;
            const Segment& previous = segments_.back();
            if (wait > duration_ - previous.start) break;
            const double next = previous.start + wait;
            if (!(next > previous.start)) {
                throw std::overflow_error("anisotropic muonium jump time cannot advance");
            }
            if (segments_.size() - 1 >= max_hops) {
                throw std::length_error("realized anisotropic muonium jump count exceeds max_hops");
            }

            rho = previous.evolution.density(next - previous.start);
            const double u = rng.Rndm() * total_rate;
            if (u < parameters.hop_rate) {
                current_field = sample_field();
            } else if (u < parameters.hop_rate + parameters.spin_exchange_rate) {
                rho = anisotropic_muonium_detail::randomize_electron_spin(rho);
            } else {
                current_axis = anisotropic_muonium_detail::sample_unit_sphere(rng);
            }

            segments_.push_back({
                next, current_field, current_axis,
                AnisotropicMuoniumEvolution(
                    rho, parameters.external_field + current_field,
                    parameters.hyperfine_iso, parameters.dipolar_D, current_axis)
            });
            if (next == duration_) break;
        }
    }

    DensityMatrix density(double time) const {
        require_nonnegative(time, "anisotropic muonium query time must be finite and >= 0");
        if (time > duration_) {
            throw std::invalid_argument("query exceeds prepared anisotropic muonium duration");
        }
        const auto after = std::upper_bound(segments_.begin(), segments_.end(), time,
            [](double t, const Segment& s) { return t < s.start; });
        const auto& segment = *(after - 1);
        return segment.evolution.density(time - segment.start);
    }

    Vec3 polarization(double time) const { return muon_polarization(density(time)); }
    const std::vector<Segment>& segments() const { return segments_; }
    std::size_t hop_count() const { return segments_.size() - 1; }
    double duration() const { return duration_; }

private:
    double duration_;
    std::vector<Segment> segments_;
};

// --- Two-state chemical exchange master equation (Cox et al., JPCM 23, 315801 (2011)) ---
// Exact 30x30 real traceless Pauli-basis Liouvillian evolution for two coupled 4x4 (mu+ (x) e-)
// radical configurations exchanging reversibly at rates W12 (1 -> 2) and W21 (2 -> 1) while
// preserving unpaired electron spin coherence and memory across transitions.
namespace two_state_detail {

using Mat30 = std::array<double, 30 * 30>;
using Vec30 = std::array<double, 30>;

inline constexpr int idx15(int alpha, int beta) {
    return 4 * alpha + beta - 1;
}

inline Mat30 mat30_identity() {
    Mat30 I{};
    for (int i = 0; i < 30; ++i) I[30 * i + i] = 1.0;
    return I;
}

inline Mat30 mat30_mul(const Mat30& A, const Mat30& B) {
    Mat30 C{};
    for (int i = 0; i < 30; ++i) {
        for (int k = 0; k < 30; ++k) {
            const double aik = A[30 * i + k];
            if (aik == 0.0) continue;
            for (int j = 0; j < 30; ++j) {
                C[30 * i + j] += aik * B[30 * k + j];
            }
        }
    }
    return C;
}

inline Vec30 mat30_vec(const Mat30& A, const Vec30& x) {
    Vec30 y{};
    for (int i = 0; i < 30; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 30; ++j) {
            sum += A[30 * i + j] * x[j];
        }
        if (!std::isfinite(sum)) {
            throw std::overflow_error("two-state exchange matrix-vector product overflow");
        }
        y[i] = sum;
    }
    return y;
}

inline double mat30_norm1(const Mat30& A) {
    double max_col = 0.0;
    for (int j = 0; j < 30; ++j) {
        double col = 0.0;
        for (int i = 0; i < 30; ++i) {
            const double v = A[30 * i + j];
            if (!std::isfinite(v)) {
                throw std::overflow_error("non-finite entry in two-state exchange matrix");
            }
            col += std::abs(v);
        }
        if (col > max_col) max_col = col;
    }
    if (!std::isfinite(max_col)) {
        throw std::overflow_error("two-state exchange matrix 1-norm overflow");
    }
    return max_col;
}

// Solves A * X = B for 30x30 real matrices using Gaussian elimination with partial pivoting.
inline Mat30 mat30_solve(Mat30 A, Mat30 B) {
    for (int k = 0; k < 30; ++k) {
        int pivot = k;
        double max_abs = std::abs(A[30 * k + k]);
        for (int i = k + 1; i < 30; ++i) {
            const double val = std::abs(A[30 * i + k]);
            if (val > max_abs) {
                max_abs = val;
                pivot = i;
            }
        }
        if (!std::isfinite(max_abs) || max_abs == 0.0) {
            throw std::overflow_error("singular matrix in two-state Pade linear solve");
        }
        if (pivot != k) {
            for (int j = k; j < 30; ++j) {
                std::swap(A[30 * k + j], A[30 * pivot + j]);
            }
            for (int j = 0; j < 30; ++j) {
                std::swap(B[30 * k + j], B[30 * pivot + j]);
            }
        }
        const double inv_diag = 1.0 / A[30 * k + k];
        for (int i = k + 1; i < 30; ++i) {
            const double factor = A[30 * i + k] * inv_diag;
            if (factor == 0.0) continue;
            A[30 * i + k] = 0.0;
            for (int j = k + 1; j < 30; ++j) {
                A[30 * i + j] -= factor * A[30 * k + j];
            }
            for (int j = 0; j < 30; ++j) {
                B[30 * i + j] -= factor * B[30 * k + j];
            }
        }
    }
    Mat30 X{};
    for (int i = 29; i >= 0; --i) {
        const double inv_diag = 1.0 / A[30 * i + i];
        for (int j = 0; j < 30; ++j) {
            double sum = B[30 * i + j];
            for (int k = i + 1; k < 30; ++k) {
                sum -= A[30 * i + k] * X[30 * k + j];
            }
            const double xij = sum * inv_diag;
            if (!std::isfinite(xij)) {
                throw std::overflow_error("non-finite solution in two-state Pade linear solve");
            }
            X[30 * i + j] = xij;
        }
    }
    return X;
}

// Higham (2005) scaling-and-squaring Pade(13, 13) matrix exponential for 30x30 real matrices.
inline Mat30 mat30_expm_pade13(const Mat30& A) {
    constexpr double b[14] = {
        64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
        1187353796428800.0,  129060195264000.0,   10559470521600.0,
        670442572800.0,      33522128640.0,       1323241920.0,
        40840800.0,          960960.0,            16380.0,
        182.0,               1.0
    };
    constexpr double theta13 = 5.371920351148152;
    const double norm1 = mat30_norm1(A);
    if (norm1 == 0.0) return mat30_identity();

    int s = 0;
    if (norm1 > theta13) {
        const double ratio = std::log2(norm1 / theta13);
        if (!std::isfinite(ratio) || ratio > 62.0) {
            throw std::overflow_error("two-state exchange matrix exponential scaling overflow");
        }
        s = std::max(0, static_cast<int>(std::ceil(ratio)));
    }
    const double inv_scale = std::ldexp(1.0, -s);
    Mat30 As{};
    for (int i = 0; i < 900; ++i) As[i] = A[i] * inv_scale;

    const Mat30 A2 = mat30_mul(As, As);
    const Mat30 A4 = mat30_mul(A2, A2);
    const Mat30 A6 = mat30_mul(A2, A4);
    const Mat30 I = mat30_identity();

    Mat30 W1{}, W2{}, Z1{}, Z2{};
    for (int i = 0; i < 900; ++i) {
        W1[i] = b[13] * A6[i] + b[11] * A4[i] + b[9] * A2[i];
        W2[i] = b[7] * A6[i] + b[5] * A4[i] + b[3] * A2[i] + b[1] * I[i];
        Z1[i] = b[12] * A6[i] + b[10] * A4[i] + b[8] * A2[i];
        Z2[i] = b[6] * A6[i] + b[4] * A4[i] + b[2] * A2[i] + b[0] * I[i];
    }
    const Mat30 A6W1 = mat30_mul(A6, W1);
    const Mat30 A6Z1 = mat30_mul(A6, Z1);
    Mat30 U_inner{}, V{};
    for (int i = 0; i < 900; ++i) {
        U_inner[i] = A6W1[i] + W2[i];
        V[i] = A6Z1[i] + Z2[i];
    }
    const Mat30 U = mat30_mul(As, U_inner);

    Mat30 denom{}, numer{};
    for (int i = 0; i < 900; ++i) {
        denom[i] = V[i] - U[i];
        numer[i] = V[i] + U[i];
    }
    Mat30 R = mat30_solve(denom, numer);
    for (int k = 0; k < s; ++k) {
        R = mat30_mul(R, R);
    }
    for (int i = 0; i < 900; ++i) {
        if (!std::isfinite(R[i])) {
            throw std::overflow_error("two-state exchange matrix exponential squaring overflow");
        }
    }
    return R;
}

// Builds the exact skew-symmetric 15x15 commutator matrix M for dc/dt = M c
// in the traceless Pauli basis B_{alpha,beta} = sigma_alpha^mu (x) sigma_beta^e.
inline std::array<double, 15 * 15> build_pauli_commutator_15(
        const Vec3& field, double A_iso, double D, const Vec3& axis_unit) {
    if (!finite(field) || !std::isfinite(field.R())) {
        throw std::invalid_argument("two-state exchange magnetic field must be finite");
    }
    if (!std::isfinite(A_iso) || !std::isfinite(D)) {
        throw std::invalid_argument("two-state exchange hyperfine A_iso and D must be finite");
    }
    if (!finite(axis_unit) || !std::isfinite(axis_unit.R()) || axis_unit.R() <= 0) {
        throw std::invalid_argument("two-state exchange symmetry axis must be finite and non-zero");
    }
    const Vec3 n = axis_unit / axis_unit.R();
    const double nv[3] = {n.X(), n.Y(), n.Z()};
    double T[3][3]{};
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            T[a][b] = ((a == b) ? (A_iso - D / 3.0) : 0.0) + D * nv[a] * nv[b];
        }
    }
    const double w_mu[3] = {-gamma_muon * field.X(), -gamma_muon * field.Y(), -gamma_muon * field.Z()};
    const double w_e[3] = {gamma_electron * field.X(), gamma_electron * field.Y(), gamma_electron * field.Z()};
    for (int a = 0; a < 3; ++a) {
        if (!std::isfinite(w_mu[a]) || !std::isfinite(w_e[a])) {
            throw std::overflow_error("two-state exchange Zeeman frequency overflow");
        }
    }

    auto levi_civita = [](int i, int j, int k) -> double {
        if ((i == 0 && j == 1 && k == 2) || (i == 1 && j == 2 && k == 0) || (i == 2 && j == 0 && k == 1)) return 1.0;
        if ((i == 0 && j == 2 && k == 1) || (i == 2 && j == 1 && k == 0) || (i == 1 && j == 0 && k == 2)) return -1.0;
        return 0.0;
    };

    std::array<double, 15 * 15> M{};
    for (int alpha = 0; alpha < 4; ++alpha) {
        for (int beta = 0; beta < 4; ++beta) {
            if (alpha == 0 && beta == 0) continue;
            const int row = idx15(alpha, beta);
            if (alpha > 0) {
                for (int k = 0; k < 3; ++k) {
                    for (int m = 0; m < 3; ++m) {
                        const double e = levi_civita(k, alpha - 1, m);
                        if (e != 0.0) M[15 * row + idx15(m + 1, beta)] -= w_mu[k] * e;
                    }
                }
            }
            if (beta > 0) {
                for (int k = 0; k < 3; ++k) {
                    for (int m = 0; m < 3; ++m) {
                        const double e = levi_civita(k, beta - 1, m);
                        if (e != 0.0) M[15 * row + idx15(alpha, m + 1)] -= w_e[k] * e;
                    }
                }
            }
            if (alpha == 0 && beta > 0) {
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        for (int m = 0; m < 3; ++m) {
                            const double e = levi_civita(b, beta - 1, m);
                            if (e != 0.0) M[15 * row + idx15(a + 1, m + 1)] -= 0.5 * T[a][b] * e;
                        }
                    }
                }
            } else if (alpha > 0 && beta == 0) {
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        for (int m = 0; m < 3; ++m) {
                            const double e = levi_civita(a, alpha - 1, m);
                            if (e != 0.0) M[15 * row + idx15(m + 1, b + 1)] -= 0.5 * T[a][b] * e;
                        }
                    }
                }
            } else {
                for (int a = 0; a < 3; ++a) {
                    for (int m = 0; m < 3; ++m) {
                        const double e = levi_civita(a, alpha - 1, m);
                        if (e != 0.0) M[15 * row + idx15(m + 1, 0)] -= 0.5 * T[a][beta - 1] * e;
                    }
                }
                for (int b = 0; b < 3; ++b) {
                    for (int m = 0; m < 3; ++m) {
                        const double e = levi_civita(b, beta - 1, m);
                        if (e != 0.0) M[15 * row + idx15(0, m + 1)] -= 0.5 * T[alpha - 1][b] * e;
                    }
                }
            }
        }
    }
    return M;
}

inline std::pair<double, std::array<double, 15>> density_to_pauli(const DensityMatrix& rho) {
    using C = std::complex<double>;
    const C pauli[4][2][2] = {
        {{C(1, 0), C(0, 0)}, {C(0, 0), C(1, 0)}},
        {{C(0, 0), C(1, 0)}, {C(1, 0), C(0, 0)}},
        {{C(0, 0), C(0, -1)}, {C(0, 1), C(0, 0)}},
        {{C(1, 0), C(0, 0)}, {C(0, 0), C(-1, 0)}},
    };
    const double trace = (rho[0][0] + rho[1][1] + rho[2][2] + rho[3][3]).real();
    std::array<double, 15> c{};
    for (int alpha = 0; alpha < 4; ++alpha) {
        for (int beta = 0; beta < 4; ++beta) {
            if (alpha == 0 && beta == 0) continue;
            C sum(0, 0);
            for (int am = 0; am < 2; ++am) {
                for (int ae = 0; ae < 2; ++ae) {
                    const int i = 2 * am + ae;
                    for (int bm = 0; bm < 2; ++bm) {
                        const C pm = pauli[alpha][bm][am];
                        if (pm == C(0, 0)) continue;
                        for (int be = 0; be < 2; ++be) {
                            const C pe = pauli[beta][be][ae];
                            if (pe == C(0, 0)) continue;
                            sum += rho[i][2 * bm + be] * pm * pe;
                        }
                    }
                }
            }
            c[idx15(alpha, beta)] = sum.real();
        }
    }
    return {trace, c};
}

inline DensityMatrix pauli_to_density(double trace, const double* c15) {
    using C = std::complex<double>;
    const C pauli[4][2][2] = {
        {{C(1, 0), C(0, 0)}, {C(0, 0), C(1, 0)}},
        {{C(0, 0), C(1, 0)}, {C(1, 0), C(0, 0)}},
        {{C(0, 0), C(0, -1)}, {C(0, 1), C(0, 0)}},
        {{C(1, 0), C(0, 0)}, {C(0, 0), C(-1, 0)}},
    };
    DensityMatrix rho{};
    for (int alpha = 0; alpha < 4; ++alpha) {
        for (int beta = 0; beta < 4; ++beta) {
            const double coeff = (alpha == 0 && beta == 0) ? trace : c15[idx15(alpha, beta)];
            if (coeff == 0.0) continue;
            const double w = 0.25 * coeff;
            for (int am = 0; am < 2; ++am) {
                for (int bm = 0; bm < 2; ++bm) {
                    const C pm = pauli[alpha][am][bm];
                    if (pm == C(0, 0)) continue;
                    for (int ae = 0; ae < 2; ++ae) {
                        for (int be = 0; be < 2; ++be) {
                            const C pe = pauli[beta][ae][be];
                            if (pe == C(0, 0)) continue;
                            rho[2 * am + ae][2 * bm + be] += w * pm * pe;
                        }
                    }
                }
            }
        }
    }
    return rho;
}

} // namespace two_state_detail

struct TwoStateExchangeParameters {
    Vec3 external_field{0, 0, 0};            // Static external magnetic field [T]
    double local_field_sigma_T = 0;          // Per-component Gaussian local-field width [T]
    double hop_rate = 0;                     // Local-field jump rate [ns^-1]
    double state1_hyperfine_iso = 0;         // State 1 signed isotropic A1 [rad/ns]
    double state1_dipolar_D = 0;             // State 1 signed axial dipolar D1 = A_par - A_perp [rad/ns]
    double state2_hyperfine_iso = 0;         // State 2 signed isotropic A2 [rad/ns] (Cox 2011: 0)
    double state2_dipolar_D = 0;             // State 2 signed axial dipolar D2 [rad/ns] (Cox 2011: 0)
    double exchange_w12 = 0;                 // Chemical exchange exit rate W12 (1 -> 2) [ns^-1]
    double exchange_w21 = 0;                 // Chemical exchange return rate W21 (2 -> 1) [ns^-1]
    double spin_exchange_rate = 0;           // Optional continuous electron spin randomization rate [ns^-1]
    double reorientation_rate = 0;           // Optional molecular axis jump rate [ns^-1]
    Vec3 fixed_axis{0, 0, 0};                // Zero => random isotropic axis per event (powder average)
    std::optional<double> initial_p1{};      // Initial State 1 fraction in [0, 1]; empty => equilibrium W21/(W12+W21)
};

class TwoStateExchangeTrajectory;

// Single-segment exact Liouvillian master-equation evolution for two exchanging radical states.
class TwoStateExchangeEvolution {
public:
    TwoStateExchangeEvolution(
            const Vec3& p0, const Vec3& field,
            double A1, double D1, double A2, double D2,
            double W12, double W21, const Vec3& axis,
            double spin_exchange_rate = 0.0,
            std::optional<double> initial_p1 = std::nullopt) {
        (void)SpinEvolution(SpinState::Muonium, p0, {});
        validate_rates(W12, W21, spin_exchange_rate, initial_p1);
        const double w_tot = W12 + W21;
        double p1_0 = 1.0;
        if (initial_p1.has_value()) {
            p1_0 = *initial_p1;
        } else if (w_tot > 0.0) {
            p1_0 = W21 / w_tot;
        }
        const double p2_0 = 1.0 - p1_0;
        std::array<double, 15> u1{}, u2{};
        u1[two_state_detail::idx15(1, 0)] = p1_0 * p0.X();
        u1[two_state_detail::idx15(2, 0)] = p1_0 * p0.Y();
        u1[two_state_detail::idx15(3, 0)] = p1_0 * p0.Z();
        u2[two_state_detail::idx15(1, 0)] = p2_0 * p0.X();
        u2[two_state_detail::idx15(2, 0)] = p2_0 * p0.Y();
        u2[two_state_detail::idx15(3, 0)] = p2_0 * p0.Z();
        init_from_pauli(p1_0, u1, p2_0, u2, field, A1, D1, A2, D2, W12, W21, axis, spin_exchange_rate);
    }

    std::pair<double, double> populations(double time) const {
        require_nonnegative(time, "two-state exchange time must be finite and >= 0");
        const double w_tot = w12_ + w21_;
        if (time == 0.0 || w_tot == 0.0) return {p1_0_, p2_0_};
        const double p_tot = p1_0_ + p2_0_;
        const double p1_eq = p_tot * (w21_ / w_tot);
        const double decay = std::exp(-w_tot * time);
        const double p1 = p1_eq + (p1_0_ - p1_eq) * decay;
        return {p1, p_tot - p1};
    }

    std::pair<DensityMatrix, DensityMatrix> state_densities(double time) const {
        const auto [p1, p2] = populations(time);
        const auto u = evolve_pauli(time);
        return {
            two_state_detail::pauli_to_density(p1, u.data()),
            two_state_detail::pauli_to_density(p2, u.data() + 15)
        };
    }

    DensityMatrix density(double time) const {
        const auto [r1, r2] = state_densities(time);
        DensityMatrix tot{};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                tot[i][j] = r1[i][j] + r2[i][j];
            }
        }
        return tot;
    }

    Vec3 polarization(double time) const {
        const auto u = evolve_pauli(time);
        constexpr int ix = two_state_detail::idx15(1, 0);
        constexpr int iy = two_state_detail::idx15(2, 0);
        constexpr int iz = two_state_detail::idx15(3, 0);
        return {u[ix] + u[15 + ix], u[iy] + u[15 + iy], u[iz] + u[15 + iz]};
    }

private:
    friend class TwoStateExchangeTrajectory;

    TwoStateExchangeEvolution(
            double p1_0, const std::array<double, 15>& u1_0,
            double p2_0, const std::array<double, 15>& u2_0,
            const Vec3& field,
            double A1, double D1, double A2, double D2,
            double W12, double W21, const Vec3& axis,
            double spin_exchange_rate) {
        validate_rates(W12, W21, spin_exchange_rate, std::nullopt);
        init_from_pauli(p1_0, u1_0, p2_0, u2_0, field, A1, D1, A2, D2, W12, W21, axis, spin_exchange_rate);
    }

    static void validate_rates(
            double W12, double W21, double spin_exchange_rate,
            std::optional<double> initial_p1) {
        require_nonnegative(W12, "two-state exchange W12 must be finite and >= 0");
        require_nonnegative(W21, "two-state exchange W21 must be finite and >= 0");
        require_nonnegative(spin_exchange_rate, "two-state spin_exchange_rate must be finite and >= 0");
        if (!std::isfinite(W12 + W21)) {
            throw std::invalid_argument("two-state exchange total rate W12 + W21 overflow");
        }
        if (initial_p1.has_value()) {
            if (!std::isfinite(*initial_p1) || *initial_p1 < 0.0 || *initial_p1 > 1.0) {
                throw std::invalid_argument("two-state initial_p1 must be finite and in [0, 1]");
            }
        }
    }

    void init_from_pauli(
            double p1_0, const std::array<double, 15>& u1_0,
            double p2_0, const std::array<double, 15>& u2_0,
            const Vec3& field,
            double A1, double D1, double A2, double D2,
            double W12, double W21, const Vec3& axis,
            double spin_exchange_rate) {
        const auto M1 = two_state_detail::build_pauli_commutator_15(field, A1, D1, axis);
        const auto M2 = two_state_detail::build_pauli_commutator_15(field, A2, D2, axis);
        w12_ = W12;
        w21_ = W21;
        p1_0_ = p1_0;
        p2_0_ = p2_0;

        const bool use_balanced = (W12 > 0.0 && W21 > 0.0);
        if (use_balanced) {
            const double w_tot = W12 + W21;
            s1_ = std::sqrt(W21 / w_tot);
            s2_ = std::sqrt(W12 / w_tot);
        } else {
            s1_ = 1.0;
            s2_ = 1.0;
        }
        const double c12 = use_balanced ? std::sqrt(W12 * W21) : W21;
        const double c21 = use_balanced ? std::sqrt(W12 * W21) : W12;

        L_.fill(0.0);
        for (int i = 0; i < 15; ++i) {
            for (int j = 0; j < 15; ++j) {
                L_[30 * i + j] = M1[15 * i + j];
                L_[30 * (15 + i) + (15 + j)] = M2[15 * i + j];
            }
            const bool has_electron_pauli = ((i + 1) % 4) != 0;
            const double se = has_electron_pauli ? spin_exchange_rate : 0.0;
            L_[30 * i + i] -= (W12 + se);
            L_[30 * (15 + i) + (15 + i)] -= (W21 + se);
            L_[30 * i + (15 + i)] = c12;
            L_[30 * (15 + i) + i] = c21;
        }
        for (int i = 0; i < 15; ++i) {
            u0_scaled_[i] = u1_0[i] / s1_;
            u0_scaled_[15 + i] = u2_0[i] / s2_;
        }
        (void)two_state_detail::mat30_norm1(L_);
    }

    two_state_detail::Vec30 evolve_pauli(double time) const {
        require_nonnegative(time, "two-state exchange time must be finite and >= 0");
        two_state_detail::Vec30 u{};
        if (time == 0.0) {
            for (int i = 0; i < 15; ++i) {
                u[i] = s1_ * u0_scaled_[i];
                u[15 + i] = s2_ * u0_scaled_[15 + i];
            }
            return u;
        }
        two_state_detail::Mat30 Lt{};
        for (int i = 0; i < 900; ++i) {
            const double v = L_[i] * time;
            if (!std::isfinite(v)) {
                throw std::overflow_error("two-state exchange L*t overflow");
            }
            Lt[i] = v;
        }
        const auto expLt = two_state_detail::mat30_expm_pade13(Lt);
        const auto u_scaled = two_state_detail::mat30_vec(expLt, u0_scaled_);
        for (int i = 0; i < 15; ++i) {
            u[i] = s1_ * u_scaled[i];
            u[15 + i] = s2_ * u_scaled[15 + i];
        }
        return u;
    }

    double w12_ = 0.0;
    double w21_ = 0.0;
    double p1_0_ = 1.0;
    double p2_0_ = 0.0;
    double s1_ = 1.0;
    double s2_ = 1.0;
    two_state_detail::Mat30 L_{};
    two_state_detail::Vec30 u0_scaled_{};
};

// Powder-averaged two-state chemical exchange trajectory with optional local-field hopping
// and molecular symmetry-axis reorientation. Chemical exchange (W12, W21) and continuous
// electron spin relaxation are integrated exactly inside each segment via exp(L * t).
class TwoStateExchangeTrajectory {
public:
    struct Segment {
        double start;
        Vec3 local_field;
        Vec3 symmetry_axis;
        TwoStateExchangeEvolution evolution;
    };

    TwoStateExchangeTrajectory(
            const Vec3& p0, const TwoStateExchangeParameters& parameters,
            double duration, TRandom3& rng, std::size_t max_hops = 100000)
        : duration_(duration) {
        // Validate all deterministic inputs before consuming any RNG draws.
        require_nonnegative(duration, "two-state exchange duration must be finite and >= 0");
        require_nonnegative(parameters.local_field_sigma_T, "local field sigma must be finite and >= 0");
        require_nonnegative(parameters.hop_rate, "hop rate must be finite and >= 0");
        require_nonnegative(parameters.reorientation_rate, "reorientation rate must be finite and >= 0");
        if (!finite(parameters.fixed_axis) || !std::isfinite(parameters.fixed_axis.R())) {
            throw std::invalid_argument("two-state exchange fixed_axis must be finite");
        }
        const bool use_fixed_axis = parameters.fixed_axis.R() > 0;
        const Vec3 probe_axis = use_fixed_axis
            ? (parameters.fixed_axis / parameters.fixed_axis.R())
            : Vec3{0, 0, 1};
        (void)TwoStateExchangeEvolution(
            p0, parameters.external_field,
            parameters.state1_hyperfine_iso, parameters.state1_dipolar_D,
            parameters.state2_hyperfine_iso, parameters.state2_dipolar_D,
            parameters.exchange_w12, parameters.exchange_w21, probe_axis,
            parameters.spin_exchange_rate, parameters.initial_p1);

        const double total_jump_rate = parameters.hop_rate + parameters.reorientation_rate;
        if (!std::isfinite(total_jump_rate)) {
            throw std::invalid_argument("two-state exchange jump rate overflow");
        }
        if (total_jump_rate * duration > static_cast<double>(max_hops)) {
            throw std::length_error("expected two-state exchange jump count exceeds max_hops");
        }

        const auto sample_field = [&]() -> Vec3 {
            if (parameters.local_field_sigma_T == 0) return {0, 0, 0};
            const double s = parameters.local_field_sigma_T;
            const Vec3 f{rng.Gaus(0, s), rng.Gaus(0, s), rng.Gaus(0, s)};
            if (!finite(f) || !std::isfinite((parameters.external_field + f).R())) {
                throw std::overflow_error("sampled two-state exchange field overflow");
            }
            return f;
        };

        Vec3 current_axis = use_fixed_axis
            ? probe_axis
            : anisotropic_muonium_detail::sample_unit_sphere(rng);
        Vec3 current_field = sample_field();
        segments_.push_back({
            0.0, current_field, current_axis,
            TwoStateExchangeEvolution(
                p0, parameters.external_field + current_field,
                parameters.state1_hyperfine_iso, parameters.state1_dipolar_D,
                parameters.state2_hyperfine_iso, parameters.state2_dipolar_D,
                parameters.exchange_w12, parameters.exchange_w21, current_axis,
                parameters.spin_exchange_rate, parameters.initial_p1)
        });

        if (total_jump_rate == 0.0 || duration == 0.0) return;

        while (true) {
            const double wait = rng.Exp(1.0) / total_jump_rate;
            const Segment& previous = segments_.back();
            if (wait > duration_ - previous.start) break;
            const double next = previous.start + wait;
            if (!(next > previous.start)) {
                throw std::overflow_error("two-state exchange jump time cannot advance");
            }
            if (segments_.size() - 1 >= max_hops) {
                throw std::length_error("realized two-state exchange jump count exceeds max_hops");
            }

            const double dt = next - previous.start;
            const auto [p1_now, p2_now] = previous.evolution.populations(dt);
            const auto u_now = previous.evolution.evolve_pauli(dt);
            std::array<double, 15> u1_now{}, u2_now{};
            for (int i = 0; i < 15; ++i) {
                u1_now[i] = u_now[i];
                u2_now[i] = u_now[15 + i];
            }

            const double r = rng.Rndm() * total_jump_rate;
            if (r < parameters.hop_rate) {
                current_field = sample_field();
            } else {
                current_axis = anisotropic_muonium_detail::sample_unit_sphere(rng);
            }

            segments_.push_back({
                next, current_field, current_axis,
                TwoStateExchangeEvolution(
                    p1_now, u1_now, p2_now, u2_now,
                    parameters.external_field + current_field,
                    parameters.state1_hyperfine_iso, parameters.state1_dipolar_D,
                    parameters.state2_hyperfine_iso, parameters.state2_dipolar_D,
                    parameters.exchange_w12, parameters.exchange_w21, current_axis,
                    parameters.spin_exchange_rate)
            });
            if (next == duration_) break;
        }
    }

    DensityMatrix density(double time) const {
        require_nonnegative(time, "two-state exchange query time must be finite and >= 0");
        if (time > duration_) {
            throw std::invalid_argument("query exceeds prepared two-state exchange duration");
        }
        const auto after = std::upper_bound(segments_.begin(), segments_.end(), time,
            [](double t, const Segment& s) { return t < s.start; });
        const auto& segment = *(after - 1);
        return segment.evolution.density(time - segment.start);
    }

    Vec3 polarization(double time) const {
        require_nonnegative(time, "two-state exchange query time must be finite and >= 0");
        if (time > duration_) {
            throw std::invalid_argument("query exceeds prepared two-state exchange duration");
        }
        const auto after = std::upper_bound(segments_.begin(), segments_.end(), time,
            [](double t, const Segment& s) { return t < s.start; });
        const auto& segment = *(after - 1);
        return segment.evolution.polarization(time - segment.start);
    }

    const std::vector<Segment>& segments() const { return segments_; }
    std::size_t hop_count() const { return segments_.size() - 1; }
    double duration() const { return duration_; }

private:
    double duration_;
    std::vector<Segment> segments_;
};

} // namespace stopped_muon
