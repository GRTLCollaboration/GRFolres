/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"
#include "CCZ4StateVariables.hpp"

#include <array>
#include <string>

// Minimal state for the modified-gravity tests: the CCZ4 variables plus the
// scalar field phi and its conjugate momentum Pi (matches the BinaryBH4dST
// example and the GRChombo FourDerivScalarTensor tests).
enum
{
    c_phi = NUM_CCZ4_VARS,
    c_Pi,

    NUM_VARS
};

namespace StateVariables
{
static constexpr int num_user_variables =
    static_cast<int>(NUM_VARS) - static_cast<int>(NUM_CCZ4_VARS);

static const amrex::Vector<std::string> additional_names = {"phi", "Pi"};

static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4StateVariables::names, additional_names);

static const std::array<BCParity, num_user_variables> user_variable_parities = {
    BCParity::even, BCParity::even};

static const std::array<BCParity, NUM_VARS> parities = ArrayTools::concatenate(
    CCZ4StateVariables::parities, user_variable_parities);

static const std::array<amrex::Real, num_user_variables>
    user_variable_asymptotic_values{};

static const std::array<amrex::Real, NUM_VARS> asymptotic_values =
    ArrayTools::concatenate(CCZ4StateVariables::asymptotic_values,
                            user_variable_asymptotic_values);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
