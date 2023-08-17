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

    c_rhoGB,

    c_JGB1,
    c_JGB2,
    c_JGB3,

    c_Weyl4_Re,
    c_Weyl4_Im,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1",     "Mom2",    "Mom3",

    "rhoGB",

    "JGB1",     "JGB2",    "JGB3",

    "Weyl4_Re", "Weyl4_Im"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */