from textwrap import dedent

import numpy as np
from pytest import mark, param
from openmc_mcnp_adapter import mcnp_str_to_model


# Cell 2 is the lattice cell (universe 50, material 8). Universes 2 through 8
# are filled with materials 1 through 7, universe 9 is void.
_TEMPLATE = dedent("""\
    lattice
    1    0          -900   fill=50
    2    8 -1.0     {region}   lat={lat} u=50
    {fill}
    10   1 -1.0     -901   u=2
    11   2 -1.0     -901   u=3
    12   3 -1.0     -901   u=4
    13   4 -1.0     -901   u=5
    14   5 -1.0     -901   u=6
    15   6 -1.0     -901   u=7
    16   7 -1.0     -901   u=8
    17   0          -901   u=9
    20   0           900

    {surfaces}
    900  so 100.0
    901  so 1000.0

    m1   1001.80c   1.0
    m2   8016.80c   1.0
    m3   26056.80c  1.0
    m4   2004.80c   1.0
    m5   3007.80c   1.0
    m6   5011.80c   1.0
    m7   13027.80c  1.0
    m8   92235.80c  1.0
    """)

# Material found in each fill universe; universe 50 is the lattice cell itself
_MATERIAL = {2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: None, 50: 8}

# Asymmetric 5x5 array with the first index varying fastest
_FILL_2D = [9, 9, 2, 3, 9,
            9, 4, 5, 6, 9,
            7, 8, 50, 2, 3,
            9, 4, 6, 9, 9,
            9, 5, 7, 9, 9]

# Three 3x3 layers, each with its own universes, with the first index varying
# fastest
_FILL_3D = [2, 3, 2, 2, 2, 2, 2, 2, 2,
            4, 4, 4, 5, 50, 4, 4, 4, 4,
            6, 6, 6, 6, 6, 6, 7, 6, 6]

_RECT_PLANES = """\
201  px  1.0
202  px -1.0
203  py  1.5
204  py -1.5
205  pz  2.0
206  pz -2.0"""

SQRT3 = 1.7320508076

# Manual Listing 10.8: hexagon with an apothem of 1 cm centered at the origin
_HEX_PLANES = """\
301  px  1.0
302  px -1.0
303  p   1.0  1.7320508076  0.0  2.0
304  p  -1.0  1.7320508076  0.0  2.0
305  p   1.0  1.7320508076  0.0 -2.0
306  p  -1.0  1.7320508076  0.0 -2.0"""


def _fill_card(indices, univ_ids, per_line=9):
    """Build the FILL card of a lattice cell as continuation lines."""
    lines = ['     fill=' + indices]
    for i in range(0, len(univ_ids), per_line):
        ids = ' '.join(str(u) for u in univ_ids[i:i + per_line])
        lines.append('     ' + ids)
    return '\n'.join(lines)


def _convert(lat, region, surfaces, indices, univ_ids):
    """Convert the template with the given lattice cell and FILL array."""
    mcnp_str = _TEMPLATE.format(
        lat=lat, region=region, surfaces=surfaces,
        fill=_fill_card(indices, univ_ids))
    return mcnp_str_to_model(mcnp_str)


def _assert_elements(model, ranges, center, vectors, univ_ids):
    """Check the material found at the center of every MCNP lattice element."""
    (i1, i2), (j1, j2), (k1, k2) = ranges
    a1, a2, a3 = (np.asarray(v, dtype=float) for v in vectors)
    n = 0
    for k in range(k1, k2 + 1):
        for j in range(j1, j2 + 1):
            for i in range(i1, i2 + 1):
                point = np.asarray(center, dtype=float) + i*a1 + j*a2 + k*a3
                cell = model.geometry.find(point)[-1]
                found = None if cell.fill is None else cell.fill.id
                assert found == _MATERIAL[univ_ids[n]], \
                    'element [{},{},{}] at {}'.format(i, j, k, point)
                n += 1


@mark.parametrize("lat,surfaces,region,indices,fill,center,vectors", [
    # The first pair of surfaces gives the direction of the first index
    param(1, _RECT_PLANES, '-203 204 -201 202', '-2:2 -2:2 0:0', _FILL_2D,
          (0., 0., 0.), ((0., 3., 0.), (2., 0., 0.), (0., 0., 0.)),
          id='rect-y-planes-first'),
    # A macrobody lists its facets in MCNP order
    param(1, '210  rpp  -1 1  -1.5 1.5  -2 2', '-210', '-2:2 -2:2 0:0',
          _FILL_2D, (0., 0., 0.), ((2., 0., 0.), (0., 3., 0.), (0., 0., 0.)),
          id='rect-rpp'),
    # The z planes listed first, bottom plane first, put the first index
    # along -z
    param(1, _RECT_PLANES, '206 -205 -201 202 -203 204', '-1:1 -1:1 -1:1',
          _FILL_3D, (0., 0., 0.), ((0., 0., -4.), (2., 0., 0.), (0., 3., 0.)),
          id='rect-3d'),
    param(2, _HEX_PLANES, '-301 302 -303 305 -304 306', '-2:2 -2:2 0:0',
          _FILL_2D, (0., 0., 0.), ((2., 0., 0.), (1., SQRT3, 0.), (0., 0., 0.)),
          id='hex-planes'),
    param(2, '1  rhp  2 0.5 -5   0 0 10   1 0 0', '-1', '-2:2 -2:2 0:0',
          _FILL_2D, (2., 0.5, 0.),
          ((2., 0., 0.), (1., SQRT3, 0.), (0., 0., 10.)),
          id='hex-rhp-off-center'),
    # Height vector of the prism pointing along -z
    param(2, '1  rhp  0 0 5   0 0 -10   1 0 0', '-1', '-1:1 -1:1 -1:1',
          _FILL_3D, (0., 0., 0.),
          ((2., 0., 0.), (1., -SQRT3, 0.), (0., 0., -10.)), id='hex-3d'),
])
def test_lattice(lat, surfaces, region, indices, fill, center, vectors):
    model = _convert(lat, region, surfaces, indices, fill)
    ranges = [tuple(map(int, pair.split(':'))) for pair in indices.split()]
    _assert_elements(model, ranges, center, vectors, fill)


@mark.parametrize("inner_first", [True, False])
def test_nested_lattice(inner_first):
    # 3x3 outer lattice (pitch 3) whose elements are universe 2, itself a 3x3
    # lattice (pitch 1) with universe 3 at its center and universe 4 elsewhere
    outer = "2 0 -2 LAT=1 U=1 FILL=-1:1 -1:1 0:0 2 2 2 2 2 2 2 2 2"
    inner = "3 0 -3 LAT=1 U=2 FILL=-1:1 -1:1 0:0 4 4 4 4 3 4 4 4 4"
    lattices = f"{inner}\n    {outer}" if inner_first else f"{outer}\n    {inner}"
    mcnp_str = dedent(f"""
    title
    1 0 -1 FILL=1
    {lattices}
    4 1 -1.0 -4 U=3
    5 2 -2.0 +4 U=3
    6 2 -2.0 -5 U=4

    1 rpp -4.5 4.5 -4.5 4.5 -4.5 4.5
    2 rpp -1.5 1.5 -1.5 1.5 -1.5 1.5
    3 rpp -0.5 0.5 -0.5 0.5 -1.5 1.5
    4 so 0.2
    5 so 10.0

    m1   1001.80c  1.0
    m2   1002.80c  1.0
    """)
    geometry = mcnp_str_to_model(mcnp_str).geometry

    # Center of the inner center element of the outer element at (3, 0, 0)
    assert geometry.find((3.0, 0.0, 0.0))[-1].fill.id == 1
    # One inner pitch over is the surrounding universe 4
    assert geometry.find((4.0, 0.0, 0.0))[-1].fill.id == 2
    # Same for the outer element at (0, 3, 0)
    assert geometry.find((0.0, 3.0, 0.0))[-1].fill.id == 1


# Lattice cell 10 cm high from z = -3 to 7 whose universe 2 holds a sphere
# offset from the origin of the universe
_LAYER_TEMPLATE = dedent("""\
    single axial layer
    1    0          -900   fill=50
    2    3 -1.0     -1   lat={lat} u=50 fill=-1:1 -1:1 {k}:{k}
         2 2 2 2 50 2 2 2 2
    10   1 -1.0     -901   u=2
    11   2 -1.0      901   u=2
    20   0           900

    {surfaces}
    900  so 100.0
    901  s 0.3 0.0 1.0 0.8

    m1   1001.80c   1.0
    m2   8016.80c   1.0
    m3   26056.80c  1.0
    """)


@mark.parametrize("lat,surfaces,k,neighbors", [
    param(1, '1  rpp  -1 1  -1.5 1.5  -3 7', 2,
          [(2., 0.), (0., 3.), (-2., -3.)], id='rect'),
    param(2, '1  rhp  0 0 -3   0 0 10   1 0 0', -3,
          [(2., 0.), (1., SQRT3), (-1., SQRT3)], id='hex'),
])
def test_single_layer_lattice(lat, surfaces, k, neighbors):
    # The single layer at axial index k becomes a 2D lattice whose universes
    # are translated to the layer; a wrong axial position of the layer would
    # miss the sphere
    mcnp_str = _LAYER_TEMPLATE.format(lat=lat, surfaces=surfaces, k=k)
    geometry = mcnp_str_to_model(mcnp_str).geometry
    assert geometry.get_all_lattices()[50].ndim == 2

    def material_at(point):
        return geometry.find(point)[-1].fill.id

    # Center of the sphere and a point above it in the elements next to the
    # [0,0,k] element, which is filled with the material of the lattice cell
    for x, y in neighbors:
        assert material_at((x + 0.3, y, 10.*k + 1.0)) == 1
        assert material_at((x + 0.3, y, 10.*k + 3.0)) == 2
    assert material_at((0.3, 0.0, 10.*k + 1.0)) == 3