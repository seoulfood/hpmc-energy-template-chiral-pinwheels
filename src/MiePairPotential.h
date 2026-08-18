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

    TODO: Rename the "MiePairPotential" to a class name that represents your potential.
    "ExamplePair" appears many times in C++, CMakeLists, and Python files. Consider using a global
    search and replace tool.
 */
class MiePairPotential : public PairPotential
    {
    public:
    MiePairPotential(std::shared_ptr<SystemDefinition> sysdef);
    virtual ~MiePairPotential() { }

    virtual LongReal energy(const LongReal r_squared,
                            const vec3<LongReal>& r_ij,
                            const unsigned int type_i,
                            const quat<LongReal>& q_i,
                            const LongReal charge_i,
                            const unsigned int type_j,
                            const quat<LongReal>& q_j,
                            const LongReal charge_j) const;

    // Compute the non-additive cuttoff radius
    // LongReal computeRCutNonAdditive(unsigned int type_i, unsigned int type_j) const;

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
            m_epsilon = 0;
            m_sigma = 0;
            m_n = 0;
            m_m = 0;
            m_r_squared = 0;

        }

        /// Construct a parameter set from a dictionary.
        ParamType(pybind11::dict params)
        {
            auto epsilon(params["epsilon"].cast<LongReal>());
            auto sigma(params["sigma"].cast<LongReal>());
            auto m(params["m"].cast<LongReal>());
            auto n(params["n"].cast<LongReal>());
            auto r_cut(params["r_cut"].cast<LongReal>());

            m_epsilon = epsilon;
            m_sigma = sigma;
            m_m = m;
            m_n = n;
            m_r_squared = r_cut * r_cut;
            m_C_epsilon = (n/(n-m)) * fast::pow((n/m), (m/(n-m))) * epsilon;

        }

        /// Convert a parameter set to a dictionary.
        pybind11::dict asDict()
        {
            pybind11::dict result;

            result["epsilon"] = m_epsilon;
            result["sigma"] = m_sigma;
            result["m"] = m_m;
            result["n"] = m_n;
            result["r_cut"] = slow::sqrt(m_r_squared);

            return result;
        }

        // TODO: rename or add parameters as needed
        LongReal m_epsilon;
        LongReal m_sigma;
        LongReal m_n;
        LongReal m_m;
        LongReal m_r_squared;

        LongReal m_C_epsilon;
        };

    /// Parameters per type pair.
    std::vector<ParamType> m_params;
    };

namespace detail
    {
//! Export the ExampleUpdater class to python
void export_MiePairPotential(pybind11::module& m);

    } // end namespace detail

    } // end namespace hpmc
    } // end namespace hoomd
