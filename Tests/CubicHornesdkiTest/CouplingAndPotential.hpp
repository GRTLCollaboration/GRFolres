/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include "simd.hpp"

class CouplingAndPotential
{
  public:
    struct params_t
    {
        double lambda_GB; // Gauss-Bonnet coupling
        double cutoff_GB; // cutoff for switching off the Gauss-Bonnet terms
                          // inside the BH
        double factor_GB; // factor for the function smoothening the GB cutoff
        double g2;        // coupling to the square of the kinetic term
        double g3;
        mutable double scalar_mass;
    };

  private:
    params_t m_params;

  public:
    template <class data_t>
    ALWAYS_INLINE data_t V(const data_t phi, const data_t X) const
    {
        return 0.5 * pow(m_params.scalar_mass * phi, 2.);
    } // V
    template <class data_t>
    ALWAYS_INLINE data_t G2(const data_t phi, const data_t X) const
    {
        return m_params.g2 * X * X;
    } // G2
    template <class data_t>
    ALWAYS_INLINE data_t dV_dphi(const data_t phi, const data_t X) const
    {
        return m_params.scalar_mass * m_params.scalar_mass * phi;
    } // dV_dphi
    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi(const data_t phi, const data_t X) const
    {
        return 0.;
    } // dG2_dphi
    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX(const data_t phi, const data_t X) const
    {
        return 2. * m_params.g2 * X;
    } // dG2_dX
    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX(const data_t phi, const data_t X) const
    {
        return 2. * m_params.g2;
    } // d2G2_dXX
    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXphi(const data_t phi, const data_t X) const
    {
        return 0.;
    } // d2G2_dXphi

    template <class data_t>
    ALWAYS_INLINE data_t G3(const data_t phi, const data_t X) const
    {
        return m_params.g3 * X;
    } // G3
    template <class data_t>
    ALWAYS_INLINE data_t dG3_dphi(const data_t phi, const data_t X) const
    {
        return 0.;
    } // dG3_dphi
    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX(const data_t phi, const data_t X) const
    {
        return m_params.g3;
    } // dG3_dX
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX(const data_t phi, const data_t X) const
    {
        return 0.;
    } // d2G3_dXX
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXphi(const data_t phi, const data_t X) const
    {
        return 0.;
    } // d2G3_dXphi
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dphiphi(const data_t phi, const data_t X) const
    {
        return 0.;
    } // d2G3_dphiphi

    //! The constructor
    CouplingAndPotential(params_t a_params) : m_params(a_params) {}

    //! Set the EsGB coupling function fhere
    template <class data_t, template <typename> class vars_t>
    void compute_coupling_and_potential(data_t &dfdphi, data_t &d2fdphi2,
                                        data_t &g2, data_t &dg2dphi,
                                        data_t &V_of_phi, data_t &dVdphi,
                                        const vars_t<data_t> &vars,
                                        const Coordinates<data_t> &coords) const
    {
        // excision setting the coupling to 0 in the interior of the BH with a
        // smooth function
        data_t cutoff_factor =
            1. + exp(-m_params.factor_GB * (vars.chi - m_params.cutoff_GB));

        // The first derivative of the GB coupling function
        dfdphi = m_params.lambda_GB / cutoff_factor;
        // The second derivative of the GB coupling function
        d2fdphi2 = M_PI * m_params.lambda_GB / cutoff_factor;
        // The coupling to the square of the kinetic term
        g2 = m_params.g2;
        // The first derivative of the g2 coupling
        dg2dphi = M_PI;
        // The potential of the scalar field
        V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0);
        // The first derivative of the potential
        dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;
    }
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
