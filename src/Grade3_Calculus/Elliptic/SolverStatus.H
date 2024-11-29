/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2024 Thomas Jefferson University and Arizona State University
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/MUON-CFD/somar.
 ******************************************************************************/
#ifndef ___EllipticSolverStatus_H__INCLUDED___
#define ___EllipticSolverStatus_H__INCLUDED___

namespace Elliptic {


/**
 * @interface SolverStatus
 * @brief     All solvers relay their convergence status to the caller
 *            via this interface.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 */
class SolverStatus
{
public:
    /**
     * @{
     * @name  Exit status constants
     * @brief I'm doing this without an enum so that child classes can add
     *        to the list if needed.
     */
    static constexpr int UNDEFINED = -1;
    static constexpr int DIVERGED  = 0;
    static constexpr int CONVERGED = 1;
    static constexpr int SINGULAR  = 2;
    static constexpr int MAXITERS  = 3;
    static constexpr int HANG      = 4;
    /** @} */

    /// Weak constructor
    SolverStatus()
    : m_statusInt(UNDEFINED)
    , m_initResNorm(-1.0)
    , m_finalResNorm(-1.0)
    {}

    /// Destructor
    virtual ~SolverStatus()
    {}

    ///
    virtual void
    clear()
    {
        m_statusInt    = UNDEFINED;
        m_initResNorm  = -1.0;
        m_finalResNorm = -1.0;
    }

    /// Sets the exit status.
    virtual void
    setSolverStatus(const int a_statusInt)
    {
        m_statusInt = a_statusInt;
    }

    /// Returns the exit status.
    virtual int
    getSolverStatus() const
    {
        return m_statusInt;
    }

    /// Convert the exit status to a readable message.
    virtual const char*
    getExitMessage() const
    {
        switch (m_statusInt) {
            case UNDEFINED:
                return "undefined";
            case DIVERGED:
                return "diverging";
            case CONVERGED:
                return "converged";
            case SINGULAR:
                return "matrix is singular";
            case MAXITERS:
                return "max iters reached";
            case HANG:
                return "hanging";
        }
        return "unknown";
    }

    /// Basic accessor
    virtual Real
    getInitResNorm() const
    {
        return m_initResNorm;
    }

    /// Basic accessor
    virtual void
    setInitResNorm(const Real& a_resNorm)
    {
        m_initResNorm = a_resNorm;
    }

    /// Basic accessor
    virtual Real
    getFinalResNorm() const
    {
        return m_finalResNorm;
    }

    /// Basic accessor
    virtual void
    setFinalResNorm(const Real& a_resNorm)
    {
        m_finalResNorm = a_resNorm;
    }

protected:
    int  m_statusInt;
    Real m_initResNorm;
    Real m_finalResNorm;
};


} // namespace Elliptic
#endif //!___EllipticSolverStatus_H__INCLUDED___
