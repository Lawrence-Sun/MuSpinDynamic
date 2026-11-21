#pragma once

#include "Math/Vector2D.h"
#include "Math/Vector3D.h"

namespace stopped_muon {

using Vec2 = ROOT::Math::XYVector;
using Vec3 = ROOT::Math::XYZVector;

// --- system of units (base units: ns, T, rad) ---
inline constexpr double ns = 1.0;                    // nanoseconds
inline constexpr double us = 1e3 * ns;               // microseconds
inline constexpr double ms = 1e3 * us;               // milliseconds
inline constexpr double GHz = 1.0 / ns;              // gigahertz
inline constexpr double MHz = 1.0 / us;              // megahertz
inline constexpr double kHz = 1.0 / ms;              // kilohertz
inline constexpr double Hz = 1e-3 * kHz;             // hertz
inline constexpr double tesla = 1.0;                 // tesla
inline constexpr double mT = 1e-3 * tesla;           // millitesla
inline constexpr double uT = 1e-3 * mT;              // microtesla
inline constexpr double nT = 1e-3 * uT;              // nanotesla
inline constexpr double rad = 1.0;                   // radians
inline constexpr double pi = 3.14159265358979323846; // pi
inline constexpr double deg = pi / 180.0;            // degrees

// --- physics constants ---
inline constexpr double muon_mass_MeV = 105.6583755;
inline constexpr double electron_mass_MeV = 0.51099895069;
inline constexpr double muon_lifetime = 2196.9811 * ns;
inline constexpr double gamma_muon = 0.850723 * rad / (ns * tesla);
inline constexpr double gamma_electron = 176.085962784 * rad / (ns * tesla);
inline constexpr double gamma_muonium = (gamma_muon + gamma_electron) / 2.0;
inline constexpr double omega_muonium_hfs = 28.04375 * rad / ns;

} // namespace stopped_muon

namespace mu_spin_dynamic = stopped_muon;
