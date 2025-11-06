# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Template HOOMD-blue component."""

import hoomd
import hoomd.conftest
from hoomd.hpmc_energy import _hpmc_energy

import inspect
from hoomd.data.typeconverter import OnlyIf, to_type_converter


class ExampleExternal(hoomd.hpmc.external.External):
    """An example external potential for HPMC.

    TODO: document the class.

    """

    _cpp_class_name = 'ExampleExternalPotential'
    _ext_module = _hpmc_energy

    def __init__(self):
        params = hoomd.data.typeparam.TypeParameter(
            'params',
            'particle_types',
            # TODO: Add entries to TypeParameterDict matching those in the C++
            # ParamType.
            hoomd.data.parameterdicts.TypeParameterDict(epsilon=float, len_keys=1),
        )
        self._add_typeparam(params)


@hoomd.logging.modify_namespace(('hpmc', 'pair', 'ExamplePair'))
class ExamplePair(hoomd.hpmc.pair.Pair):
    r"""An example pair potential for HPMC.

    Args:
        default_r_cut (float): Default cutoff radius :math:`[\mathrm{length}]`.

    TODO: document the class.

    """

    _cpp_class_name = 'ExamplePairPotential'
    _ext_module = _hpmc_energy

    def __init__(self, default_r_cut=None):
        if default_r_cut is None:
            default_r_cut = float
        else:
            default_r_cut = float(default_r_cut)

        params = hoomd.data.typeparam.TypeParameter(
            'params',
            'particle_types',
            # TODO: Add entries to TypeParameterDict matching those in the C++
            # ParamType.
            hoomd.data.parameterdicts.TypeParameterDict(
                A=float, B=float, r_cut=default_r_cut, len_keys=2
            ),
        )
        self._add_typeparam(params)


@hoomd.logging.modify_namespace(('hpmc', 'pair', 'HarmonicAnglePair'))
class HarmonicAnglePair(hoomd.hpmc.pair.Pair):
    r"""An example pair potential for HPMC.

    Args:
        default_r_cut (float): Default cutoff radius :math:`[\mathrm{length}]`.

    TODO: document the class.

    """

    _cpp_class_name = 'HarmonicAnglePairPotential'
    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(hoomd.hpmc.pair.Pair._doc_inherited)
    )
    _ext_module = _hpmc_energy

    def __init__(self, default_theta = None, default_equivalent_theta=None, default_k = None, default_r_cut = None):
        if default_r_cut is None:
            default_r_cut = float
        else:
            default_r_cut = float(default_r_cut)

        if default_theta is None:
            default_theta = float
        else:
            default_theta = float(default_theta)

        if default_equivalent_theta is None:
            default_equivalent_theta = float
        else:
            default_equivalent_theta = float(default_equivalent_theta)

        if default_k is None:
            default_k = float
        else:
            default_k = float(default_k)

        params = hoomd.data.typeparam.TypeParameter(
            "params",
            "particle_types",
            hoomd.data.parameterdicts.TypeParameterDict(
                theta=default_theta,
                equivalent_theta = default_equivalent_theta,
                k=default_k,
                r_cut=default_r_cut,
                len_keys=2,
            ),
        )
        self._add_typeparam(params)


