/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CCZ4RHS.hpp"
#include "Coordinates.hpp"
#include "CouplingAndPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourDerivScalarTensor.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedGauge.hpp"
#include "Tensor.hpp"
#include <iostream>

template <class data_t> struct Vars : public CCZ4::Vars<data_t>
{
    data_t phi;  // the scalar field
    data_t Pi; // the curvature of the scalar field

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        CCZ4::Vars<data_t>::enum_mapping(mapping_function);
        VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
        VarsTools::define_enum_mapping(mapping_function, c_Pi, Pi);
    }
};

typedef ModifiedCCZ4RHS<FourDerivScalarTensor<CouplingAndPotential>, MovingPunctureGauge,
                     FourthOrderDerivatives, ModifiedGauge> MyModifiedGravityClass;

int main(int argc, char *argv[])
{
    int failed = 0;

    CCZ4::params_t m_params;
    ModifiedGauge::params_t m_params_modified_gauge;
    CouplingAndPotential::params_t m_params_coupling_and_potential;    

    MyModifiedGravityClass::Vars<double> rhs;
    MyModifiedGravityClass::Vars<double> vars;
    MyModifiedGravityClass::Vars<Tensor<1, double>> d1;
    MyModifiedGravityClass::Diff2Vars<Tensor<2, double>> d2;
    MyModifiedGravityClass::Vars<double> advec;

#include "values1.hpp" //Including the auto generated file with values    
    
    double dx = 1.;
    double sigma = 1.;
    double cosmological_constant = 0.;

    CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> obj(
        m_params, dx, sigma, 0, cosmological_constant);
    obj.rhs_equation(rhs, vars, d1, d2, advec);

    const IndexTM<int, 3> ind = {0, 0, 0};
    IntVect vect(ind);
    Coordinates<double> coords(vect, dx, {0., 0., 0.});
    
    ModifiedGauge modified_gauge(m_params_modified_gauge);
    CouplingAndPotential coupling_and_potential(m_params_coupling_and_potential);
    FourDerivScalarTensor<CouplingAndPotential> fdst(coupling_and_potential);
    MyModifiedGravityClass my_modified_ccz4(fdst, m_params, modified_gauge, dx, sigma, 0, 1./(16.*M_PI));

    // add functions a(x) and b(x) of the modified gauge
    my_modified_ccz4.add_a_b_rhs<double>(rhs, vars, d1, d2, advec, coords);

    // add RHS matter terms from EM Tensor
    my_modified_ccz4.add_emtensor_rhs<double>(rhs, vars, d1, d2, advec,
                                               coords);

    // add evolution of matter fields themselves
    fdst.add_matter_rhs<double>(rhs, vars, d1, d2, advec);

    // solve linear system for the matter fields that require it (e.g. EsGB)
    fdst.solve_lhs<double>(rhs, vars, d1, d2, advec);

    // Compare
    double diff;
    diff = rhs.lapse - dlapsedt_known;
    if (diff > 1e-10 or diff < -1e-10)
    {
        std::cout << "RHS of lapse wrong" << std::endl;
        std::cout << "value: " << rhs.lapse << std::endl;
        std::cout << "correct value: " << dlapsedt_known << std::endl;
        std::cout << "diff: " << diff << std::endl;
        failed = -1;
    }
    FOR(i)
    {
        diff = rhs.shift[i] - dshiftdt_known[i];
        if (diff > 1e-10 or diff < -1e-10)
        {
            std::cout << "RHS of shift wrong in component [" << i << "]"
                      << std::endl;
            std::cout << "value: " << rhs.shift[i] << std::endl;
            std::cout << "correct value: " << dshiftdt_known[i] << std::endl;
            std::cout << "diff: " << diff << std::endl;
            failed = -1;
        }
    }
    FOR(i)
    {
        diff = rhs.B[i] - dBdt_known[i];
        if (diff > 1e-10 or diff < -1e-10)
        {
            std::cout << "RHS of B wrong in component [" << i << "]"
                      << std::endl;
            std::cout << "value: " << rhs.B[i] << std::endl;
            std::cout << "correct value: " << dBdt_known[i] << std::endl;
            std::cout << "diff: " << diff << std::endl;
            failed = -1;
        }
    }
    FOR(i, j)
    {
        diff = rhs.h[i][j] - dhdt_known[i][j];
        if (diff > 1e-10 or diff < -1e-10)
        {
            std::cout << "RHS of h wrong in component [" << i << "][" << j
                      << "]" << std::endl;
            std::cout << "value: " << rhs.h[i][j] << std::endl;
            std::cout << "correct value: " << dhdt_known[i][j] << std::endl;
            std::cout << "diff: " << diff << std::endl;
            failed = -1;
        }
    }
    diff = rhs.chi - dchidt_known;
    if (diff > 1e-10 or diff < -1e-10)
    {
        std::cout << "RHS of chi wrong" << std::endl;
        std::cout << "value: " << rhs.chi << std::endl;
        std::cout << "correct value: " << dchidt_known << std::endl;
        std::cout << "diff: " << diff << std::endl;
        failed = -1;
    }
    diff = rhs.phi - dphidt_known;
    if (diff > 1e-10 or diff < -1e-10)
    {
        std::cout << "RHS of phi wrong" << std::endl;
        std::cout << "value: " << rhs.phi << std::endl;
        std::cout << "correct value: " << dphidt_known << std::endl;
        std::cout << "diff: " << diff << std::endl;
        failed = -1;
    }
    diff = rhs.Theta - dThetadt_known;
    if (diff > 1e-10 or diff < -1e-10)
    {
        std::cout << "RHS of Theta wrong" << std::endl;
        std::cout << "value: " << rhs.Theta << std::endl;
        std::cout << "correct value: " << dThetadt_known << std::endl;
        std::cout << "diff: " << diff << std::endl;
        failed = -1;
    }
    FOR(i)
    {
        diff = rhs.Gamma[i] - dGammadt_known[i];
        if (diff > 1e-10 or diff < -1e-10)
        {
            std::cout << "RHS of Gamma wrong in component [" << i << "]"
                      << std::endl;
            std::cout << "value: " << rhs.Gamma[i] << std::endl;
            std::cout << "correct value: " << dGammadt_known[i] << std::endl;
            std::cout << "diff: " << diff << std::endl;
            failed = -1;
        }
    }
    diff = rhs.K - dKdt_known;
    if (diff > 1e-10 or diff < -1e-10)
    {
        std::cout << "RHS of K wrong" << std::endl;
        std::cout << "value: " << rhs.K << std::endl;
        std::cout << "correct value: " << dKdt_known << std::endl;
        std::cout << "diff: " << diff << std::endl;
        failed = -1;
    }
    diff = rhs.Pi - dPidt_known;
    if (diff > 1e-10 or diff < -1e-10)
    {
        std::cout << "RHS of Pi wrong" << std::endl;
        std::cout << "value: " << rhs.Pi << std::endl;
        std::cout << "correct value: " << dPidt_known << std::endl;
        std::cout << "diff: " << diff << std::endl;
        failed = -1;
    }
    FOR(i, j)
    {
        diff = rhs.A[i][j] - dAdt_known[i][j];
        if (diff > 1e-10 or diff < -1e-10)
        {
            std::cout << "RHS of A wrong in component [" << i << "][" << j
                      << "]" << std::endl;
            std::cout << "value: " << rhs.A[i][j] << std::endl;
            std::cout << "correct value: " << dAdt_known[i][j] << std::endl;
            std::cout << "diff: " << diff << std::endl;
            failed = -1;
        }
    }
    if (failed == 0)
        std::cout << "FourDerivScalarTensor test passed..." << std::endl;
    else
        std::cout << "FourDerivScalarTensor test failed..." << std::endl;

    return failed;
}
