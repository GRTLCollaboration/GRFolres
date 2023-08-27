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
        double V_of_phi;
        double dVdphi;
        double G2;
        double dG2dphi;
        double dG2dX;
        double d2G2dXX;
        double d2G2dXphi;
        double G3;
        double dG3dphi;
        double dG3dX;
        double d2G3dXX;
        double d2G3dXphi;
        double d2G3dphiphi;
    };

  private:
    params_t m_params;

  public:
    template <class data_t>
    ALWAYS_INLINE data_t V(const data_t phi, const data_t X) const
    {
        return m_params.V_of_phi;
    } // V
    template <class data_t>
    ALWAYS_INLINE data_t G2(const data_t phi, const data_t X) const
    {
        return m_params.G2;
    } // G2
    template <class data_t>
    ALWAYS_INLINE data_t dV_dphi(const data_t phi, const data_t X) const
    {
        return m_params.dVdphi;
    } // dV_dphi
    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi(const data_t phi, const data_t X) const
    {
        return m_params.dG2dphi;
    } // dG2_dphi
    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX(const data_t phi, const data_t X) const
    {
        return m_params.dG2dX;
    } // dG2_dX
    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX(const data_t phi, const data_t X) const
    {
        return m_params.d2G2dXX;
    } // d2G2_dXX
    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXphi(const data_t phi, const data_t X) const
    {
        return m_params.d2G2dXphi;
    } // d2G2_dXphi

    template <class data_t>
    ALWAYS_INLINE data_t G3(const data_t phi, const data_t X) const
    {
        return m_params.G3;
    } // G3
    template <class data_t>
    ALWAYS_INLINE data_t dG3_dphi(const data_t phi, const data_t X) const
    {
        return m_params.dG3dphi;
    } // dG3_dphi
    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX(const data_t phi, const data_t X) const
    {
        return m_params.dG3dX;
    } // dG3_dX
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX(const data_t phi, const data_t X) const
    {
        return m_params.d2G3dXX;
    } // d2G3_dXX
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXphi(const data_t phi, const data_t X) const
    {
        return m_params.d2G3dXphi;
    } // d2G3_dXphi
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dphiphi(const data_t phi, const data_t X) const
    {
        return m_params.d2G3dphiphi;
    } // d2G3_dphiphi

    //! The constructor
    CouplingAndPotential(params_t a_params) : m_params(a_params) {}
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
