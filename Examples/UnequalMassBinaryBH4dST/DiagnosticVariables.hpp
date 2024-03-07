/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_rho_phi,
    c_rho_g2,
    c_rho_g3,
    c_rho_GB,

    c_weak_coupling_condition_g2,
    c_weak_coupling_condition_g3,
    c_weak_coupling_condition_GB,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1",
    "Mom2",
    "Mom3",

    "Weyl4_Re",
    "Weyl4_Im",

    "rho_phi",
    "rho_g2",
    "rho_g3",
    "rho_GB",

    "weak_coupling_condition_g2",
    "weak_coupling_condition_g3",
    "weak_coupling_condition_GB"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
