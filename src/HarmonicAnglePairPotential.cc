// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "HarmonicAnglePairPotential.h"

namespace hoomd
    {
namespace hpmc
    {

HarmonicAnglePairPotential::HarmonicAnglePairPotential(std::shared_ptr<SystemDefinition> sysdef)
    : PairPotential(sysdef), m_params(m_type_param_index.getNumElements())
    {
    }

LongReal HarmonicAnglePairPotential::energy(const LongReal r_squared,
                                      const vec3<LongReal>& r_ij,
                                      const unsigned int type_i,
                                      const quat<LongReal>& q_i,
                                      const LongReal charge_i,
                                      const unsigned int type_j,
                                      const quat<LongReal>& q_j,
                                      const LongReal charge_j) const
    {
    unsigned int param_index = m_type_param_index(type_i, type_j);
    const auto& param = m_params[param_index];

    vec3<LongReal> v1_i = rotate(q_i, vec3<LongReal>(1, 0, 0));
    vec3<LongReal> v2_i = rotate(q_i, vec3<LongReal>(0, 1, 0));

    vec3<LongReal> v1_j = rotate(q_j, vec3<LongReal>(1, 0, 0));
    vec3<LongReal> v2_j = rotate(q_j, vec3<LongReal>(0, 1, 0));

    Scalar theta_1_ij = acos(dot(v1_i, v1_j));
    Scalar r_ij_mag = dot(r_ij, r_ij);

    vec3<LongReal> proj_rij_v2_i = v2_i - ((dot(v2_i, r_ij)/r_ij_mag)*r_ij);
    vec3<LongReal> proj_rij_v2_j = v2_j - ((dot(v2_j, r_ij)/r_ij_mag)*r_ij);

    Scalar theta_2_ij = acos(dot(proj_rij_v2_i, proj_rij_v2_j));

    Scalar theta_theta_star_squared = fast::pow(theta_2_ij-param.m_theta_star, 2);

    //Scalar invr_rsq = 1 / r_squared;
    LongReal energy = - 0.5 * param.m_k * (theta_theta_star_squared);
    return energy;

    }

void HarmonicAnglePairPotential::setParamsPython(pybind11::tuple particle_types, pybind11::dict params)
    {
    auto pdata = m_sysdef->getParticleData();
    auto type_i = pdata->getTypeByName(particle_types[0].cast<std::string>());
    auto type_j = pdata->getTypeByName(particle_types[1].cast<std::string>());
    unsigned int param_index_1 = m_type_param_index(type_i, type_j);
    m_params[param_index_1] = ParamType(params);
    unsigned int param_index_2 = m_type_param_index(type_j, type_i);
    m_params[param_index_2] = ParamType(params);

    notifyRCutChanged();
    }

pybind11::dict HarmonicAnglePairPotential::getParamsPython(pybind11::tuple particle_types)
    {
    auto pdata = m_sysdef->getParticleData();
    auto type_i = pdata->getTypeByName(particle_types[0].cast<std::string>());
    auto type_j = pdata->getTypeByName(particle_types[1].cast<std::string>());
    unsigned int param_index = m_type_param_index(type_i, type_j);
    return m_params[param_index].asDict();
    }

// namespace detail
//     {
// //! Export the ExampleUpdater class to python
// void export_HarmonicAnglePairPotential(pybind11::module& m)
//     {
//     pybind11::class_<HarmonicAnglePairPotential, PairPotential, std::shared_ptr<HarmonicAnglePairPotential>>(
//         m,
//         "HarmonicAnglePairPotential")
//         .def(pybind11::init<std::shared_ptr<SystemDefinition>>())
//         .def("setParams", &HarmonicAnglePairPotential::setParamsPython)
//         .def("getParams", &HarmonicAnglePairPotential::getParamsPython);

//     }

namespace detail
    {
void export_HarmonicAnglePairPotential(pybind11::module& m)
    {
    pybind11::class_<HarmonicAnglePairPotential,
                     PairPotential,
                     std::shared_ptr<HarmonicAnglePairPotential>>(m, "HarmonicAnglePairPotential")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>>())
        .def("setParams", &HarmonicAnglePairPotential::setParamsPython)
        .def("getParams", &HarmonicAnglePairPotential::getParamsPython);
    }
 

    } // end namespace detail

    } // end namespace hpmc
} // end namespace hoomd

