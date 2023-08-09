/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef KERRBH4DSTLEVEL_HPP_
#define KERRBH4DSTLEVEL_HPP_

#include "BHAMR.hpp"
#include "CouplingAndPotential.hpp"
#include "DefaultLevelFactory.hpp"
#include "FourDerivScalarTensor.hpp"
#include "GRAMRLevel.hpp"
#include "ModifiedPunctureGauge.hpp"

class KerrBH4dSTLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<KerrBH4dSTLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);

    // Typedef for 4dST
    typedef FourDerivScalarTensor<CouplingAndPotential>
        FourDerivScalarTensorWithCouplingAndPotential;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

#ifdef CH_USE_HDF5
    /// Things to do before writing a plot file
    virtual void prePlotLevel() override;
#endif /* CH_USE_HDF5 */

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells() override;

    virtual void
    computeTaggingCriterion(FArrayBox &tagging_criterion,
                            const FArrayBox &current_state) override;
};

#endif /* KERRBH4DSTLEVEL_HPP_ */
