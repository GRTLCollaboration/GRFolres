/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBHTESTFIELD4DSTLEVEL_HPP_
#define BINARYBHTESTFIELD4DSTLEVEL_HPP_

#include "BHAMR.hpp"
#include "CouplingAndPotential.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "ModifiedPunctureGauge.hpp"
#include "TestField4dST.hpp"

class BinaryBHTestField4dSTLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<BinaryBHTestField4dSTLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);

    // Typedef for 4dST
    typedef TestField4dST<CouplingAndPotential>
        TestField4dSTWithCouplingAndPotential;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells() override;

    /// Identify and tag the cells that need higher resolution
    virtual void
    computeTaggingCriterion(FArrayBox &tagging_criterion,
                            const FArrayBox &current_state) override;

    // to do post each time step on every level
    virtual void specificPostTimeStep() override;

#ifdef CH_USE_HDF5
    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;
#endif /* CH_USE_HDF5 */
};

#endif /* BINARYBHTESTFIELD4DSTLEVEL_HPP_ */
