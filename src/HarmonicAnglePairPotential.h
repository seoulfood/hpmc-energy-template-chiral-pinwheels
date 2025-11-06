// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#include <hoomd/VectorMath.h>
#include <hoomd/hpmc/PairPotential.h>
#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace hpmc
    {

/** Example pair potential for use with HPMC simulations.

    TODO: Rename the "HarmonicAnglePairPotential" to a class name that represents your potential.
    "HarmonicAnglePair" appears many times in C++, CMakeLists, and Python files. Consider using a global
    search and replace tool.
 */

class HarmonicAnglePairPotential : public hpmc::PairPotential
    {
    public:
    HarmonicAnglePairPotential(std::shared_ptr<SystemDefinition> sysdef);
    virtual ~HarmonicAnglePairPotential() { }

    virtual LongReal energy(const LongReal r_squared,
                            const vec3<LongReal>& r_ij,
                            const unsigned int type_i,
                            const quat<LongReal>& q_i,
                            const LongReal charge_i,
                            const unsigned int type_j,
                            const quat<LongReal>& q_j,
                            const LongReal charge_j) const;

    /// Compute the non-additive cuttoff radius
    virtual LongReal computeRCutNonAdditive(unsigned int type_i, unsigned int type_j) const
        {
        unsigned int param_index = m_type_param_index(type_i, type_j);
        return slow::sqrt(m_params[param_index].m_r_squared);
        }


    /// Set type-pair-dependent parameters to the potential.
    void setParamsPython(pybind11::tuple particle_types, pybind11::dict params);

    /// Get type-pair-dependent parameters.
    pybind11::dict getParamsPython(pybind11::tuple particle_types);

    protected:
    /// per-type-pair parameters
    struct ParamType
        {
        ParamType() 
        { 
            m_theta_star = 0;
            m_equiv_theta = 0;
            m_k = 0;
            m_r_squared = 0;
        }

        /// Construct a parameter set from a dictionary.
        ParamType(pybind11::dict params)
        {
            auto theta_star(params["theta"].cast<LongReal>());
            auto k(params["k"].cast<LongReal>());
            auto equiv_theta(params["equivalent_theta"].cast<LongReal>());
            auto r_cut(params["r_cut"].cast<LongReal>());
            
            m_theta_star = theta_star;
            m_equiv_theta = equiv_theta;
            m_k = k;
            m_r_squared = r_cut * r_cut;


        }

        pybind11::dict asDict()
            {
            pybind11::dict result;

            result["theta"] = m_theta_star;
            result["k"] = m_k;
            result["equivalent_theta"] = m_equiv_theta;
            result["r_cut"] = slow::sqrt(m_r_squared);

            return result;
            }

        // TODO: rename or add parameters as needed
        LongReal m_theta_star;
        LongReal m_equiv_theta;
        LongReal m_k;
        LongReal m_r_squared;
        };

    /// Parameters per type pair.
    std::vector<ParamType> m_params;
    };

namespace detail
    {
//! Export the ExampleUpdater class to python
void export_HarmonicAnglePairPotential(pybind11::module& m);

    } // end namespace detail


    } // end namespace hpmc
    } // end namespace hoomd
