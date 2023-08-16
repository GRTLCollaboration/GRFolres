
#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

class CouplingAndPotential
{
  public:
    struct params_t
    {
        // 'mutable' - mass can be changed even if object is constant
        // this is useful to set it to 0 (temporarily)
        // to calculate the WFC
        mutable double scalar_mass;
        double g3;
        double g2;
    };

    params_t m_params;

    template <class data_t>
    ALWAYS_INLINE data_t G2(const data_t phi, const data_t X) const
    {
        return m_params.g2 * X * X - 0.5 * pow(m_params.scalar_mass * phi, 2.);
    } // G2
    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi(const data_t phi, const data_t X) const
    {
        return -m_params.scalar_mass * m_params.scalar_mass * phi;
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
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
