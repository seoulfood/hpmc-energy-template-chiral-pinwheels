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
    Scalar r_ij_mag = fast::sqrt(dot(r_ij, r_ij));

    vec3<LongReal> proj_rij_v2_i = v2_i - ((dot(v2_i, r_ij)/r_ij_mag)*r_ij);
    vec3<LongReal> proj_rij_v2_j = v2_j - ((dot(v2_j, r_ij)/r_ij_mag)*r_ij);

    //std::cout << proj_rij_v2_i.x << ", " << proj_rij_v2_i.y << ", " << proj_rij_v2_i.z << std::endl;
    //std::cout << proj_rij_v2_j.x << ", " << proj_rij_v2_j.y << ", " << proj_rij_v2_j.z << std::endl;

    Scalar theta_2_ij = acos(dot(proj_rij_v2_i, proj_rij_v2_j));
    if(theta_2_ij != theta_2_ij){
        if(dot(proj_rij_v2_i, proj_rij_v2_j) > 1){
            theta_2_ij = 0.0;
        }
        else{
            theta_2_ij = M_PI;
        }
    }

    //it's calculating the raw angle correctly just fine
    //but it's adding in a chip at the third peak

    Scalar floor_output = 0;
    Scalar modulo_output = 0;

    //std::cout << "theta_2_ij: " << theta_2_ij << std::endl;
    if(param.m_equiv_theta != 0){
        theta_2_ij = 2 * abs( (theta_2_ij/param.m_equiv_theta) - (floor((theta_2_ij/param.m_equiv_theta) + (0.5)))); //copied from python??
        floor_output = floor((theta_2_ij/param.m_equiv_theta) + (0.5));
        modulo_output = fmod(theta_2_ij-(param.m_equiv_theta/4), (param.m_equiv_theta));

    }

    Scalar theta_theta_star_squared = fast::pow(theta_2_ij-param.m_theta_star, 2);

    //Scalar invr_rsq = 1 / r_squared;
    LongReal energy = - 0.5 * param.m_k * (theta_theta_star_squared);// * fast::pow(invr_rsq, 2);
    return energy;
    //return theta_theta_star_squared;
    //return theta_2_ij;
    //return param.m_equiv_theta;
    //return floor_output;
    //return modulo_output;

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

