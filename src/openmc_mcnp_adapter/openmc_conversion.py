# SPDX-FileCopyrightText: 2022-2025 UChicago Argonne, LLC and contributors
# SPDX-License-Identifier: MIT

import argparse
from math import pi, isclose, sqrt
import os
import re
import tempfile
import warnings

import numpy as np
import openmc
from openmc.data import get_thermal_name
from openmc.data.ace import get_metadata
from openmc.model.surface_composite import (
    CompositeSurface,
    RightCircularCylinder as RCC,
    RectangularParallelepiped as RPP,
    OrthogonalBox as BOX,
    ConicalFrustum as TRC,
    HexagonalPrism,
)
from openmc.model import surface_composite

from .parse import parse, _COMPLEMENT_RE, _CELL_FILL_RE


# The facet number corresponding to the SurfaceComposite's surface by
# attribute name and whether or not to flip the sense of that surface
# based on the facet surface's relationship to the composite surface region
_MACROBODY_FACETS = {
    BOX: {
        1: ('ax1_max', False),
        2: ('ax1_min', True),
        3: ('ax2_max', False),
        4: ('ax2_min', True),
        5: ('ax3_max', False),
        6: ('ax3_min', True),
    },
    RCC: {
        1: ('cyl', False),
        2: ('top', False),
        3: ('bottom', True)
    },
    HexagonalPrism: {
        1: ('plane_max', False),
        2: ('plane_min', True),
        3: ('upper_right', False),
        4: ('lower_left', True),
        5: ('upper_left', False),
        6: ('lower_right', True),
        7: ('top', False),
        8: ('bottom', True),
    },
    RPP: {
        1: ('xmax', False),
        2: ('xmin', True),
        3: ('ymax', False),
        4: ('ymin', True),
        5: ('zmax', False),
        6: ('zmin', True)
    },
    TRC: {
        1: ('cone', False),
        2: ('plane_top', False),
        3: ('plane_bottom', True),
    }
}


def rotation_matrix(v1, v2):
    """Compute rotation matrix that would rotate v1 into v2.

    Parameters
    ----------
    v1 : numpy.ndarray
        Unrotated vector
    v2 : numpy.ndarray
        Rotated vector

    Returns
    -------
    3x3 rotation matrix

    """
    # Normalize vectors and compute cosine
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)
    cos_angle = float(np.clip(np.dot(u1, u2), -1.0, 1.0))

    I = np.identity(3)

    # Handle special case where vectors are parallel or anti-parallel
    if isclose(abs(cos_angle), 1.0, rel_tol=1e-8):
        if cos_angle > 0.0:
            return I
        else:
            # Proper 180° rotation: rotate about any axis that is orthogonal
            # with u1. Because |k| = 1 and cos(180°) = -1, the rotation matrix
            # is simply K = I + 2 * (k k^T - I) = 2 k k^T - I

            # Choose reference vector not parallel to u1
            ref = np.array([1.0, 0.0, 0.0]) if abs(u1[0]) < 0.9 else np.array([0.0, 1.0, 0.0])

            # Create orthogonal unit vector
            k = np.cross(u1, ref)
            k /= np.linalg.norm(k)

            # Create rotation matrix
            return 2.0 * np.outer(k, k) - I
    else:
        # Calculate rotation angle
        sin_angle = np.sqrt(1 - cos_angle*cos_angle)

        # Calculate axis of rotation
        axis = np.cross(u1, u2)
        axis /= np.linalg.norm(axis)

        # Create cross-product matrix K
        kx, ky, kz = axis
        K = np.array([
            [0.0, -kz, ky],
            [kz, 0.0, -kx],
            [-ky, kx, 0.0]
        ])

        # Create rotation matrix using Rodrigues' rotation formula
        return I + K * sin_angle + (K @ K) * (1 - cos_angle)


def get_openmc_materials(materials, expand_elements: bool = True):
    """Get OpenMC materials from MCNP materials

    Parameters
    ----------
    materials : list
        List of MCNP material information

    Returns
    -------
    dict
        Dictionary mapping material ID to :class:`openmc.Material`

    """
    openmc_materials = {}
    for m in materials.values():
        if 'id' not in m:
            continue
        material = openmc.Material(m['id'])
        for nuclide, percent in m['nuclides']:
            if '.' in nuclide:
                zaid, xs = nuclide.split('.')
            else:
                zaid = nuclide
            name, element, Z, A, metastable = get_metadata(int(zaid), 'mcnp')
            if percent < 0:
                if (A > 0) or (not expand_elements):
                    material.add_nuclide(name, abs(percent), 'wo')
                else:
                    material.add_element(element, abs(percent), 'wo')
            else:
                if (A > 0) or (not expand_elements):
                    material.add_nuclide(name, percent, 'ao')
                else:
                    material.add_element(element, percent, 'ao')

        if 'sab' in m:
            for sab in m['sab']:
                if '.' in sab:
                    name, xs = sab.split('.')
                else:
                    name = sab
                material.add_s_alpha_beta(get_thermal_name(name))
        openmc_materials[m['id']] = material

    return openmc_materials


def get_openmc_surfaces(surfaces, data):
    """Get OpenMC surfaces from MCNP surfaces

    Parameters
    ----------
    surfaces : list
        List of MCNP surfaces
    data : dict
        MCNP data-block information

    Returns
    -------
    dict
        Dictionary mapping surface ID to :class:`openmc.Surface` instance

    """
    # Ensure that autogenerated IDs for surfaces don't conflict
    openmc.Surface.next_id = max(s['id'] for s in surfaces) + 1

    openmc_surfaces = {}
    for s in surfaces:
        coeffs = s['coefficients']
        if s['mnemonic'] == 'p':
            if len(coeffs) == 9:
                p1 = coeffs[:3]
                p2 = coeffs[3:6]
                p3 = coeffs[6:]
                surf = openmc.Plane.from_points(p1, p2, p3, surface_id=s['id'])

                # Helper function to flip signs on plane coefficients
                def flip_sense(surf):
                    surf.a = -surf.a
                    surf.b = -surf.b
                    surf.c = -surf.c
                    surf.d = -surf.d

                # Enforce MCNP sense requirements
                if surf.d != 0.0:
                    if surf.d < 0.0:
                        flip_sense(surf)
                elif surf.c != 0.0:
                    if surf.c < 0.0:
                        flip_sense(surf)
                elif surf.b != 0.0:
                    if surf.b < 0.0:
                        flip_sense(surf)
                elif surf.a != 0.0:
                    if surf.a < 0.0:
                        flip_sense(surf)
                else:
                    raise ValueError(f"Plane {s['id']} appears to be a line? ({coeffs})")
            else:
                A, B, C, D = coeffs
                surf = openmc.Plane(surface_id=s['id'], a=A, b=B, c=C, d=D)
        elif s['mnemonic'] == 'px':
            surf = openmc.XPlane(surface_id=s['id'], x0=coeffs[0])
        elif s['mnemonic'] == 'py':
            surf = openmc.YPlane(surface_id=s['id'], y0=coeffs[0])
        elif s['mnemonic'] == 'pz':
            surf = openmc.ZPlane(surface_id=s['id'], z0=coeffs[0])
        elif s['mnemonic'] == 'so':
            surf = openmc.Sphere(surface_id=s['id'], r=coeffs[0])
        elif s['mnemonic'] in ('s', 'sph'):
            x0, y0, z0, R = coeffs
            surf = openmc.Sphere(surface_id=s['id'], x0=x0, y0=y0, z0=z0, r=R)
        elif s['mnemonic'] == 'sx':
            x0, R = coeffs
            surf = openmc.Sphere(surface_id=s['id'], x0=x0, r=R)
        elif s['mnemonic'] == 'sy':
            y0, R = coeffs
            surf = openmc.Sphere(surface_id=s['id'], y0=y0, r=R)
        elif s['mnemonic'] == 'sz':
            z0, R = coeffs
            surf = openmc.Sphere(surface_id=s['id'], z0=z0, r=R)
        elif s['mnemonic'] == 'c/x':
            y0, z0, R = coeffs
            surf = openmc.XCylinder(surface_id=s['id'], y0=y0, z0=z0, r=R)
        elif s['mnemonic'] == 'c/y':
            x0, z0, R = coeffs
            surf = openmc.YCylinder(surface_id=s['id'], x0=x0, z0=z0, r=R)
        elif s['mnemonic'] == 'c/z':
            x0, y0, R = coeffs
            surf = openmc.ZCylinder(surface_id=s['id'], x0=x0, y0=y0, r=R)
        elif s['mnemonic'] == 'cx':
            surf = openmc.XCylinder(surface_id=s['id'], r=coeffs[0])
        elif s['mnemonic'] == 'cy':
            surf = openmc.YCylinder(surface_id=s['id'], r=coeffs[0])
        elif s['mnemonic'] == 'cz':
            surf = openmc.ZCylinder(surface_id=s['id'], r=coeffs[0])
        elif s['mnemonic'] in ('k/x', 'k/y', 'k/z'):
            x0, y0, z0, R2 = coeffs[:4]
            if len(coeffs) > 4 and coeffs[4] != 0.0:
                up = (coeffs[4] > 0.0)
                if s['mnemonic'] == 'k/x':
                    surf = surface_composite.XConeOneSided(x0=x0, y0=y0, z0=z0, r2=R2, up=up)
                elif s['mnemonic'] == 'k/y':
                    surf = surface_composite.YConeOneSided(x0=x0, y0=y0, z0=z0, r2=R2, up=up)
                else:
                    surf = surface_composite.ZConeOneSided(x0=x0, y0=y0, z0=z0, r2=R2, up=up)
            else:
                if s['mnemonic'] == 'k/x':
                    surf = openmc.XCone(surface_id=s['id'], x0=x0, y0=y0, z0=z0, r2=R2)
                elif s['mnemonic'] == 'k/y':
                    surf = openmc.YCone(surface_id=s['id'], x0=x0, y0=y0, z0=z0, r2=R2)
                else:
                    surf = openmc.ZCone(surface_id=s['id'], x0=x0, y0=y0, z0=z0, r2=R2)
        elif s['mnemonic'] in ('kx', 'ky', 'kz'):
            x, R2 = coeffs[:2]
            if len(coeffs) > 2 and coeffs[2] != 0.0:
                up = (coeffs[2] > 0.0)
                if s['mnemonic'] == 'kx':
                    surf = surface_composite.XConeOneSided(x0=x, r2=R2, up=up)
                elif s['mnemonic'] == 'ky':
                    surf = surface_composite.YConeOneSided(y0=x, r2=R2, up=up)
                else:
                    surf = surface_composite.ZConeOneSided(z0=x, r2=R2, up=up)
            else:
                if s['mnemonic'] == 'kx':
                    surf = openmc.XCone(surface_id=s['id'], x0=x, r2=R2)
                elif s['mnemonic'] == 'ky':
                    surf = openmc.YCone(surface_id=s['id'], y0=x, r2=R2)
                else:
                    surf = openmc.ZCone(surface_id=s['id'], z0=x, r2=R2)
        elif s['mnemonic'] == 'sq':
            a, b, c, D, E, F, G, x, y, z = coeffs
            d = e = f = 0.0
            g = 2*(D - a*x)
            h = 2*(E - b*y)
            j = 2*(F - c*z)
            k = a*x*x + b*y*y + c*z*z + 2*(D*x + E*y + F*z) + G
            surf = openmc.Quadric(surface_id=s['id'], a=a, b=b, c=c, d=d, e=e,
                                  f=f, g=g, h=h, j=j, k=k)
        elif s['mnemonic'] == 'gq':
            a, b, c, d, e, f, g, h, j, k = coeffs
            surf = openmc.Quadric(surface_id=s['id'], a=a, b=b, c=c, d=d, e=e,
                                  f=f, g=g, h=h, j=j, k=k)
        elif s['mnemonic'] in ('tx', 'ty', 'tz'):
            x0, y0, z0, a, b, c = coeffs
            if isclose(a, 0.0, abs_tol=1e-12) and isclose(b, c):
                warnings.warn(
                    f"Degenerate torus surface {s['id']} (A=0, B=C) converted "
                    f"to an openmc.Sphere of radius {b}."
                )
                surf = openmc.Sphere(surface_id=s['id'], x0=x0, y0=y0, z0=z0, r=b)
            else:
                cls = getattr(openmc, f"{s['mnemonic'][1].upper()}Torus")
                surf = cls(surface_id=s['id'], x0=x0, y0=y0, z0=z0, a=a, b=b, c=c)
        elif s['mnemonic'] in ('x', 'y', 'z'):
            axis = s['mnemonic'].upper()
            cls_plane = getattr(openmc, f'{axis}Plane')
            cls_cylinder = getattr(openmc, f'{axis}Cylinder')
            cls_cone = getattr(surface_composite, f'{axis}ConeOneSided')
            if len(coeffs) == 2:
                x1, r1 = coeffs
                surf = cls_plane(x1, surface_id=s['id'])
            elif len(coeffs) == 4:
                x1, r1, x2, r2 = coeffs
                if x1 == x2:
                    surf = cls_plane(x1, surface_id=s['id'])
                elif r1 == r2:
                    surf = cls_cylinder(r=r1, surface_id=s['id'])
                else:
                    dr = r2 - r1
                    dx = x2 - x1
                    grad = dx/dr
                    offset = x2 - grad*r2
                    angle = (-1/grad)**2

                    # decide if we want the up or down part of the
                    # cone since one sheet is used
                    up = grad >= 0
                    kwargs = {f"{s['mnemonic']}0": offset, "r2": angle, "up": up}
                    surf = cls_cone(**kwargs)
            else:
                raise NotImplementedError(f"{s['mnemonic']} surface with {len(coeffs)} parameters")
        elif s['mnemonic'] == 'rcc':
            vx, vy, vz, hx, hy, hz, r = coeffs
            if hx == 0.0 and hy == 0.0 and hz > 0.0:
                surf = RCC((vx, vy, vz), hz, r, axis='z')
            elif hy == 0.0 and hz == 0.0 and hx > 0.0:
                surf = RCC((vx, vy, vz), hx, r, axis='x')
            elif hx == 0.0 and hz == 0.0 and hy > 0.0:
                surf = RCC((vx, vy, vz), hy, r, axis='y')
            else:
                # Create vectors for Z-axis and cylinder orientation
                u = np.array([0., 0., 1.])
                h = np.array([hx, hy, hz])

                # Determine rotation matrix to transform u -> h
                rotation = rotation_matrix(u, h)

                # Create RCC aligned with Z-axis
                height = np.linalg.norm(h)
                surf = RCC((vx, vy, vz), height, r, axis='z')

                # Rotate the RCC
                surf = surf.rotate(rotation, pivot=(vx, vy, vz))

        elif s['mnemonic'] == 'rpp':
            surf = RPP(*coeffs)
        elif s['mnemonic'] == 'box':
            v = coeffs[:3]
            a1 = coeffs[3:6]
            a2 = coeffs[6:9]
            if len(coeffs) == 12:
                a3 = coeffs[9:]
                surf = BOX(v, a1, a2, a3)
            else:
                surf = BOX(v, a1, a2)
        elif s['mnemonic'] == 'trc':
            v = coeffs[:3]
            h = coeffs[3:6]
            r1 = coeffs[6]
            r2 = coeffs[7]
            surf = TRC(v, h, r1, r2)
        elif s['mnemonic'] in ('rhp', 'hex'):
            # Missing entries are zero, so that a lone value after the height
            # vector is the x component of r
            coeffs = list(coeffs)
            if len(coeffs) < 9:
                coeffs += [0.0]*(9 - len(coeffs))
            v, h, r = (np.array(coeffs[i:i + 3]) for i in (0, 3, 6))
            if len(coeffs) > 9:
                # The facet vectors s and t are only supported for a regular
                # hexagon, where they are r rotated by 60 and 120 degrees
                # about h
                coeffs += [0.0]*(15 - len(coeffs))
                axis = h/np.linalg.norm(h)
                x = r - axis*np.dot(axis, r)
                y = np.cross(axis, x)
                tolerance = 1e-5*np.linalg.norm(x)
                for vec, angle in ((coeffs[9:12], pi/3),
                                   (coeffs[12:15], 2*pi/3)):
                    regular = x*np.cos(angle) + y*np.sin(angle)
                    if not np.allclose(vec, regular, atol=tolerance):
                        raise NotImplementedError(
                            f"{s['mnemonic'].upper()} surface {s['id']} is "
                            "not a regular hexagonal prism")

            height = np.linalg.norm(h)
            if height == 0.0:
                raise ValueError(f"Height vector of {s['mnemonic'].upper()} "
                                 f"surface {s['id']} must be nonzero")
            axis = h/height

            # Only the component of r perpendicular to the axis is meaningful
            r = r - axis*np.dot(axis, r)
            apothem = np.linalg.norm(r)
            if apothem == 0.0:
                raise ValueError(f"Facet vector of {s['mnemonic'].upper()} "
                                 f"surface {s['id']} must be nonzero")

            # Regular hexagon with facets perpendicular to x and an apothem
            # equal to the length of r; the prism is infinite along its axis
            # if the height is at least 1e6 cm, as in MCNP
            ends = ({} if height >= 1.0e6
                    else {'zmin': 0.0, 'zmax': height})
            surf = HexagonalPrism(edge_length=2*apothem/sqrt(3),
                                  orientation='y', **ends)

            # Rotate x onto r and z onto the axis, then move the bottom to v
            x_axis = r/apothem
            rotation = np.column_stack((x_axis, np.cross(axis, x_axis), axis))
            surf = surf.rotate(rotation).translate(v)
        else:
            raise NotImplementedError('Surface type "{}" not supported'
                                      .format(s['mnemonic']))

        # Set boundary conditions
        boundary = s.get('boundary')
        if boundary == 'reflective':
            surf.boundary_type = 'reflective'
        elif boundary == 'white':
            surf.boundary_type = 'white'
        elif boundary == 'periodic':
            surf.boundary_type = 'periodic'

        if 'tr' in s:
            tr_num = s['tr']
            displacement, rotation = data['tr'][tr_num]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", openmc.IDWarning)
                surf = surf.translate(displacement, inplace=True)
                if rotation is not None:
                    surf = surf.rotate(rotation, pivot=displacement, inplace=True)

        openmc_surfaces[s['id']] = surf

        # For macrobodies, we also need to add generated surfaces to dictionary
        if isinstance(surf, surface_composite.CompositeSurface):
            openmc_surfaces.update((-surf).get_surfaces())

    # Make another pass to set periodic surfaces
    for s in surfaces:
        periodic_surface_id = s.get('periodic_surface')
        if periodic_surface_id is not None:
            surf.periodic_surface = openmc_surfaces[periodic_surface_id]

    return openmc_surfaces


def replace_macrobody_facets(region: str, surfaces: dict) -> str:
    """Replace macrobody facet identifiers with the corresponding OpenMC surface

    Parameters
    ----------
    region : str
        Boolean expression relating surface half-spaces.
    surfaces : dict
        Dictionary mapping surface ID to :class:`openmc.Surface`

    Returns
    -------
    str
        An updated expression replacing the macrobody facet specification with
        the ID of that surface in the SurfaceComposite object.
    """
    # Get list of facets, sorted by string length to ensure that, e.g.,
    # replacing '3.1' will not happen before replacing '23.1'
    facets = set(re.findall(r'[-+]?\d+\.\d', region))
    facets = sorted(facets, key=len, reverse=True)

    for facet in facets:
        # Break up macrobody facet into surface ID and facet number
        surface_id, facet_num = facet.split('.')
        surface_id = int(surface_id)
        facet_num = int(facet_num)

        # Get corresponding composite surface
        composite_surf = surfaces[abs(surface_id)]
        if isinstance(composite_surf, CompositeSurface):
            # Get composite surface and whether to flip sense
            facet_attr, flip_sense = _MACROBODY_FACETS[type(composite_surf)][facet_num]
            facet_surface = getattr(composite_surf, facet_attr)
        else:
            warnings.warn(f'Macrobody facet {facet} ignored (not a macrobody)')
            facet_surface = composite_surf
            flip_sense = False

        # Get corresponding surface and its ID
        facet_id = facet_surface.id

        # starting with a positive facet ID, adjust for:
        # a) the specified sense in the original expression
        if surface_id < 0:
            facet_id = -facet_id

        # b) the sense of the facet with respect to the macrobody
        if flip_sense:
            facet_id = -facet_id

        # ensure this surface is present in the surfaces dictionary
        surfaces[facet_surface.id] = facet_surface

        # re-build the region expression with this entry in place of the macrobody facet entry
        region = region.replace(facet, str(facet_id))

    return region


def _macrobody_facets(region: str, surfaces: dict) -> str:
    """Expand a region that is a single negative macrobody into its facets

    A lattice cell bounded by a macrobody takes its index directions from the
    facets in MCNP order: [1,0,0] lies beyond facet 1, [-1,0,0] beyond facet 2,
    and so on, so the facets must be listed in that order.

    Parameters
    ----------
    region : str
        Boolean expression relating surface half-spaces.
    surfaces : dict
        Dictionary mapping surface ID to :class:`openmc.Surface`

    Returns
    -------
    str
        The expression with a lone negative macrobody replaced by the
        intersection of its facets in MCNP order, otherwise unchanged.
    """
    match = re.fullmatch(r'\s*-(\d+)\s*', region)
    if match is None:
        return region
    surface_id = match.group(1)
    surf = surfaces.get(int(surface_id))
    if not isinstance(surf, CompositeSurface):
        return region
    facets = _MACROBODY_FACETS[type(surf)]
    return ' '.join(f'-{surface_id}.{num}'
                    for num, (attr, _) in sorted(facets.items())
                    if hasattr(surf, attr))


# Tolerance used when comparing lattice vectors
_LATTICE_TOL = 1.0e-6


def _new_lattice(parameters, uid):
    """Create an empty lattice matching the LAT card of an MCNP cell

    Parameters
    ----------
    parameters : dict
        Parameters of the MCNP lattice cell
    uid : int
        Universe ID of the lattice

    Returns
    -------
    openmc.Lattice
        Rectangular lattice for LAT=1, hexagonal lattice for LAT=2

    """
    if int(parameters['lat']) == 2:
        return openmc.HexLattice(uid)
    else:
        return openmc.RectLattice(uid)


def _lattice_element(region):
    """Determine the geometry of a single lattice element

    The region of a lattice cell is an intersection of planar half-spaces
    listed in pairs of opposite facets. In MCNP, the element [1,0,0] lies
    beyond the first surface listed, [-1,0,0] beyond the second, [0,1,0]
    beyond the third, and so on. Each pair of facets thus gives the
    displacement to the neighboring element across it as well as one equation
    for the center of the element.

    Parameters
    ----------
    region : openmc.Region
        Region of the lattice cell

    Returns
    -------
    center : numpy.ndarray
        Center of the [0,0,0] lattice element
    vectors : list of numpy.ndarray
        Displacement to the neighboring element across each pair of facets, in
        the order in which the facets are listed

    """
    n = len(region) if isinstance(region, openmc.Intersection) else 1
    if n < 4:
        raise NotImplementedError('One-dimensional lattices not supported')
    if n % 2 or n > 8 or \
            not all(isinstance(node, openmc.Halfspace) for node in region):
        raise ValueError('Lattice cell must be bounded by four, six or eight '
                         'planar surfaces')

    normals, offsets = [], []
    for node in region:
        surface = node.surface
        if not isinstance(surface, openmc.PlaneMixin):
            raise ValueError('Lattice cell must be bounded by planar surfaces '
                             'but surface {} is a {}'.format(
                                 surface.id, type(surface).__name__))
        normal = np.array([surface.a, surface.b, surface.c])
        length = np.linalg.norm(normal)

        # A cell on the negative side of a plane has an outward normal that
        # points along the normal of the plane
        sense = 1.0 if node.side == '-' else -1.0
        normals.append(sense*normal/length)
        offsets.append(sense*surface.d/length)

    vectors, rows, rhs = [], [], []
    for k in range(0, n, 2):
        u1, s1 = normals[k], offsets[k]
        u2, s2 = normals[k + 1], offsets[k + 1]
        if not np.allclose(u1, -u2, atol=_LATTICE_TOL):
            raise ValueError('Facets {} and {} of a lattice cell are not '
                             'parallel and opposite'.format(k + 1, k + 2))
        vectors.append((s1 + s2)*u1)
        rows.append(u1)
        rhs.append((s1 - s2)/2)

    center = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    return center, vectors


def _fill_rect_lattice(lattice, center, vectors, ranges, univ_ids,
                       get_universe):
    """Set the geometry and the universes of a rectangular lattice

    Parameters
    ----------
    lattice : openmc.RectLattice
        Lattice to be filled
    center : numpy.ndarray
        Center of the [0,0,0] lattice element
    vectors : list of numpy.ndarray
        Displacement to the neighboring element across each pair of facets, as
        given by :func:`_lattice_element`
    ranges : iterable of tuple of int
        Lower and upper index of the FILL array along each of the three MCNP
        lattice directions
    univ_ids : numpy.ndarray
        Universe IDs of the FILL array with the first index varying fastest
    get_universe : callable
        Function returning the universe with a given ID

    """
    ndim = len(vectors)
    if ndim > 3:
        raise ValueError('Rectangular lattice cell must be bounded by four or '
                         'six planar surfaces')

    # Coordinate axis along which each lattice direction lies
    axes = []
    for vec in vectors:
        axis = int(np.argmax(np.abs(vec)))
        if not np.isclose(abs(vec[axis]), np.linalg.norm(vec),
                          atol=_LATTICE_TOL):
            raise NotImplementedError('Rectangular lattices with sides not '
                                      'perpendicular to a coordinate axis are '
                                      'not supported')
        axes.append(axis)
    if ndim == 2 and sorted(axes) != [0, 1]:
        raise NotImplementedError('2D lattice with basis other than x-y not '
                                  'supported')
    if sorted(axes) != list(range(ndim)):
        raise ValueError('Lattice directions must lie along different axes')

    # Universe IDs as an array indexed by ([k], j, i)
    counts = [upper - lower + 1 for lower, upper in ranges[:ndim]]
    univ_ids = np.asarray(univ_ids).reshape(counts[::-1])

    # Make each index increase along its axis and order the array as
    # ([z], y, x)
    for m, vec in enumerate(vectors):
        if vec[axes[m]] < 0.0:
            univ_ids = np.flip(univ_ids, axis=ndim - 1 - m)
    univ_ids = np.transpose(univ_ids, [ndim - 1 - axes.index(axis)
                                       for axis in reversed(range(ndim))])

    pitch, lower_left, dimension = [], [], []
    for axis in range(ndim):
        m = axes.index(axis)
        step = vectors[m][axis]
        lower, upper = ranges[m]
        pitch.append(abs(step))
        lower_left.append(center[axis] + min(lower*step, upper*step)
                          - abs(step)/2)
        dimension.append(upper - lower + 1)

    lattice.pitch = pitch
    lattice.lower_left = lower_left
    lattice.dimension = dimension

    # Fill universes in OpenMC lattice, reversing y direction
    lattice.universes = np.vectorize(get_universe)(univ_ids)[..., ::-1, :]


def _fill_hex_lattice(lattice, center, vectors, ranges, univ_ids, get_universe,
                      filler):
    """Set the geometry and the universes of a hexagonal lattice

    Parameters
    ----------
    lattice : openmc.HexLattice
        Lattice to be filled
    center : numpy.ndarray
        Center of the [0,0,0] lattice element
    vectors : list of numpy.ndarray
        Displacement to the neighboring element across each pair of facets, as
        given by :func:`_lattice_element`
    ranges : iterable of tuple of int
        Lower and upper index of the FILL array along each of the three MCNP
        lattice directions
    univ_ids : numpy.ndarray
        Universe IDs of the FILL array with the first index varying fastest
    get_universe : callable
        Function returning the universe with a given ID
    filler : openmc.Universe
        Universe used for the hexagons that are not part of the FILL array

    """
    if len(vectors) not in (3, 4):
        raise ValueError('Hexagonal lattice cell must be bounded by six or '
                         'eight planar surfaces')
    a1, a2, a5 = vectors[:3]

    # In MCNP, [1,0,0] is across the first facet, [0,1,0] across the third and
    # [-1,1,0] across the fifth
    if not np.allclose(a5, a2 - a1, atol=_LATTICE_TOL):
        raise ValueError('Hexagonal lattice surface order does not follow the '
                         'MCNP hexagonal lattice convention (the element '
                         'across the fifth facet must be [-1,1,0])')
    if not (np.isclose(a1[2], 0.0, atol=_LATTICE_TOL) and
            np.isclose(a2[2], 0.0, atol=_LATTICE_TOL)):
        raise NotImplementedError('Only hexagonal lattices along the z-axis '
                                  'are supported')

    # Check that the hexagon is regular
    pitch = np.linalg.norm(a1)
    cos_angle = np.dot(a1, a2)/(pitch*np.linalg.norm(a2))
    if not (np.isclose(np.linalg.norm(a2), pitch, atol=_LATTICE_TOL) and
            np.isclose(abs(cos_angle), 0.5, atol=_LATTICE_TOL)):
        raise ValueError('Irregular hexagonal lattices are not supported')

    # Orientation 'x' has facets perpendicular to x, 'y' perpendicular to y
    angle = np.degrees(np.arctan2(a1[1], a1[0])) % 60.0
    if np.isclose(angle, 0.0, atol=_LATTICE_TOL) or \
            np.isclose(angle, 60.0, atol=_LATTICE_TOL):
        lattice.orientation = 'x'
        basis = np.array([[pitch, pitch/2], [0., pitch*sqrt(3)/2]])
    elif np.isclose(angle, 30.0, atol=_LATTICE_TOL):
        lattice.orientation = 'y'
        basis = np.array([[pitch*sqrt(3)/2, 0.], [pitch/2, pitch]])
    else:
        raise NotImplementedError('Rotated hexagonal lattices not supported')

    # Matrix converting MCNP indices into OpenMC (x, alpha) indices
    matrix = np.linalg.solve(basis, np.array([a1[:2], a2[:2]]).T)
    matrix = np.rint(matrix).astype(int)

    (i1, i2), (j1, j2), (k1, k2) = ranges

    # Number of rings needed to cover the four corners of the MCNP array
    corners = [matrix @ (i, j) for i in (i1, i2) for j in (j1, j2)]
    num_rings = 1 + max(max(abs(x), abs(a), abs(x + a)) for x, a in corners)

    def empty_rings():
        return [[filler]*max(6*(num_rings - 1 - r), 1) for r in range(num_rings)]

    three_d = (len(vectors) == 4)
    if three_d:
        a3 = vectors[3]
        if not np.allclose(a3[:2], 0.0, atol=_LATTICE_TOL):
            raise NotImplementedError('Only hexagonal lattices along the '
                                      'z-axis are supported')
        lattice.universes = [empty_rings() for _ in range(k2 - k1 + 1)]
        lattice.pitch = (pitch, abs(a3[2]))
        lattice.center = (center[0], center[1],
                          center[2] + (k1 + k2)/2*a3[2])
    else:
        a3 = np.zeros(3)
        lattice.universes = empty_rings()
        lattice.pitch = (pitch,)
        lattice.center = (center[0], center[1])

    n = 0
    for k in range(k1, k2 + 1):
        for j in range(j1, j2 + 1):
            for i in range(i1, i2 + 1):
                x, alpha = matrix @ (i, j)
                if three_d:
                    z = k - k1 if a3[2] > 0 else k2 - k
                    iz, ring, pos = lattice.get_universe_index((x, alpha, z))
                    lattice.universes[iz][ring][pos] = \
                        get_universe(univ_ids[n])
                else:
                    ring, pos = lattice.get_universe_index((x, alpha))
                    lattice.universes[ring][pos] = get_universe(univ_ids[n])
                n += 1


def get_openmc_universes(cells, surfaces, materials, data):
    """Get OpenMC surfaces from MCNP surfaces

    Parameters
    ----------
    cells : list
        List of MCNP cells
    surfaces : dict
        Dictionary mapping surface ID to :class:`openmc.Surface`
    materials : dict
        Dictionary mapping material ID to :class:`openmc.Material`
    data : dict
        MCNP data-block information

    Returns
    -------
    dict
        Dictionary mapping universe ID to :class:`openmc.Universe` instance

    """
    openmc_cells = {}
    cell_by_id = {c['id']: c for c in cells}
    universes = {}
    root_universe = openmc.Universe(0)
    universes[0] = root_universe

    # Determine maximum IDs so that autogenerated IDs don't conflict
    openmc.Cell.next_id = max(c['id'] for c in cells) + 1
    all_univ_ids = set()
    for c in cells:
        if 'u' in c['parameters']:
            all_univ_ids.add(abs(int(c['parameters']['u'])))
    if all_univ_ids:
        openmc.Universe.next_id = max(all_univ_ids) + 1

    # Cell-complements pose a unique challenge for conversion because the
    # referenced cell may have a region that was translated, so we can't simply
    # replace the cell-complement by what appears on the referenced
    # cell. Instead, we loop over all the cells and construct regions for all
    # cells without cell complements. Then, we handle the remaining cells by
    # replacing the cell-complement with the string representation of the actual
    # region that was already converted
    has_cell_complement = []
    translate_memo = {}
    for c in cells:
        # Skip cells that have cell-complements to be handled later
        match = _COMPLEMENT_RE.search(c['region'])
        if match:
            has_cell_complement.append(c)
            continue

        # Assign region to cell based on expression
        region = c['region'].replace('#', '~').replace(':', '|')

        # A lattice cell bounded by a macrobody needs its facets in MCNP order
        if 'lat' in c['parameters']:
            region = _macrobody_facets(region, surfaces)

        # Replace macrobody facet specifiers in the region expression
        if '.' in region:
            region = replace_macrobody_facets(region, surfaces)

        try:
            c['_region'] = openmc.Region.from_expression(region, surfaces)
        except Exception:
            raise ValueError('Could not parse region for cell (ID={}): {}'
                             .format(c['id'], region))

        if 'trcl' in c['parameters'] or '*trcl' in c['parameters']:
            if 'trcl' in c['parameters']:
                trcl = c['parameters']['trcl'].strip()
                use_degrees = False
            else:
                trcl = c['parameters']['*trcl'].strip()
                use_degrees = True

            # Apply transformation to fill
            if 'fill' in c['parameters']:
                # TODO: Check for existing transformations on the fill
                fill = c['parameters']['fill']
                if use_degrees:
                    c['parameters']['*fill'] = f'{fill} {trcl}'
                    c['parameters'].pop('fill')
                else:
                    c['parameters']['fill'] = f'{fill} {trcl}'

            if not trcl.startswith('('):
                raise NotImplementedError(
                    'TRn card not supported (cell {}).'.format(c['id']))

            # Drop parentheses
            trcl = trcl[1:-1].split()

            # Get displacement vector
            vector = np.array([float(c) for c in trcl[:3]])

            if len(trcl) > 3:
                # If displacement vector origin is -1, reverse displacement vector
                if len(trcl) == 13:
                    if int(trcl[12]) == -1:
                        vector *= -1
                c['_region'] = c['_region'].translate(vector, translate_memo)

                rotation_matrix = np.array([float(x) for x in trcl[3:12]]).reshape((3, 3))
                if use_degrees:
                    rotation_matrix = np.cos(rotation_matrix * pi/180.0)
                c['_region'] = c['_region'].rotate(rotation_matrix.T, pivot=vector)
            else:
                c['_region'] = c['_region'].translate(vector, translate_memo)

            # Update surfaces dictionary with new surfaces
            for surf_id, surf in c['_region'].get_surfaces().items():
                surfaces[surf_id] = surf
                if isinstance(surf, surface_composite.CompositeSurface):
                    surfaces.update((-surf).get_surfaces())

    has_cell_complement_ordered = []
    def add_to_ordered(c):
        region = c['region']
        matches = _COMPLEMENT_RE.findall(region)
        for _, other_id in matches:
            other_cell = cell_by_id[int(other_id)]
            if other_cell in has_cell_complement:
                add_to_ordered(other_cell)
        if c not in has_cell_complement_ordered:
            has_cell_complement_ordered.append(c)
    for c in has_cell_complement:
        add_to_ordered(c)

    # Now that all cells without cell-complements have been handled, we loop
    # over the remaining ones and convert any cell-complement expressions by
    # using str(region)
    for c in has_cell_complement_ordered:
        # Replace cell-complement with regular complement
        region = c['region']
        matches = _COMPLEMENT_RE.findall(region)
        assert matches
        for _, other_id in matches:
            other_cell = cell_by_id[int(other_id)]
            try:
                r = ~other_cell['_region']
            except KeyError:
                raise NotImplementedError(
                    'Cannot handle nested cell-complements for cell {}: {}'
                    .format(c['id'], c['region']))
            region = _COMPLEMENT_RE.sub(str(r), region, count=1)

        # Assign region to cell based on expression
        region = region.replace('#', '~').replace(':', '|')

        # A lattice cell bounded by a macrobody needs its facets in MCNP order
        if 'lat' in c['parameters']:
            region = _macrobody_facets(region, surfaces)

        # Replace macrobody facet specifiers in the region expression
        if '.' in region:
            region = replace_macrobody_facets(region, surfaces)

        try:
            c['_region'] = openmc.Region.from_expression(region, surfaces)
        except Exception:
            raise ValueError('Could not parse region for cell (ID={}): {}'
                             .format(c['id'], region))

        # assume these cells are not translated themselves
        assert 'trcl' not in c['parameters']

    # Now that all cell regions have been converted, the next loop is to create
    # actual Cell/Universe/Lattice objects
    material_clones = {}
    # Parameters of lattice cells by universe ID, along with a cache of
    # wrapper universes for them (OpenMC can't place a lattice directly inside
    # another lattice)
    lattice_params = {abs(int(ci['parameters']['u'])): ci['parameters']
                      for ci in cells
                      if 'lat' in ci['parameters'] and 'u' in ci['parameters']}
    lattice_wrappers = {}

    for c in cells:
        cell = openmc.Cell(cell_id=c['id'])

        # Assign region to cell based on expression
        cell.region = c['_region']

        # Add cell to universes if necessary
        if 'u' in c['parameters']:
            if 'lat' not in c['parameters']:
                # Note: a negative universe indicates that the cell is not
                # truncated by the boundary of a higher level cell.
                uid = abs(int(c['parameters']['u']))
                if uid not in universes:
                    universes[uid] = openmc.Universe(uid)
                universes[uid].add_cell(cell)
        else:
            root_universe.add_cell(cell)

        # Look for vacuum boundary condition
        if isinstance(cell.region, openmc.Union):
            if all([isinstance(n, openmc.Halfspace) for n in cell.region]):
                if 'imp:n' in c['parameters'] and float(c['parameters']['imp:n']) == 0.0:
                    for n in cell.region:
                        if n.surface.boundary_type == 'transmission':
                            n.surface.boundary_type = 'vacuum'
                    root_universe.remove_cell(cell)
        elif isinstance(cell.region, openmc.Halfspace):
            if 'imp:n' in c['parameters'] and float(c['parameters']['imp:n']) == 0.0:
                if cell.region.surface.boundary_type == 'transmission':
                    cell.region.surface.boundary_type = 'vacuum'
                root_universe.remove_cell(cell)

        # Determine material fill if present -- this is not assigned until later
        # in case it's used in a lattice (need to create an extra universe then)
        cell_material_id: int = c['material']
        if cell_material_id > 0:
            mat = materials[cell_material_id]
            cell_density = c['density']
            if mat.density is None:
                if cell_density > 0:
                    mat.set_density('atom/b-cm', cell_density)
                else:
                    mat.set_density('g/cm3', abs(cell_density))
            elif mat.density != abs(c['density']):
                key = (cell_material_id, cell_density)
                if key not in material_clones:
                    material_clones[key] = mat = mat.clone()
                    if c['density'] > 0:
                        mat.set_density('atom/b-cm', c['density'])
                    else:
                        mat.set_density('g/cm3', abs(c['density']))
                else:
                    mat = material_clones[key]

        # Create lattices
        if 'fill' in c['parameters'] or '*fill' in c['parameters']:
            if 'lat' in c['parameters']:
                # Cell filled with Lattice
                uid = abs(int(c['parameters']['u']))
                if uid not in universes:
                    universes[uid] = _new_lattice(c['parameters'], uid)
                lattice = universes[uid]
                hexagonal = isinstance(lattice, openmc.HexLattice)

                def get_universe(uid):
                    if uid not in universes:
                        if uid in lattice_params:
                            universes[uid] = _new_lattice(
                                lattice_params[uid], uid)
                        else:
                            universes[uid] = openmc.Universe(uid)
                    univ = universes[uid]
                    if isinstance(univ, openmc.Lattice):
                        # OpenMC cannot place a lattice directly inside
                        # another lattice, so return a universe wrapping it
                        if uid not in lattice_wrappers:
                            lattice_wrappers[uid] = openmc.Universe(
                                cells=[openmc.Cell(fill=univ)])
                        return lattice_wrappers[uid]
                    return univ

                # Geometry of a single lattice element
                center, vectors = _lattice_element(cell.region)

                # Get extent of lattice
                fill = c['parameters']['fill']
                words = fill.split()

                # If there's only a single parameter, the lattice is infinite
                inf_lattice = (len(words) == 1)

                if inf_lattice:
                    ranges = [(0, 0), (0, 0), (0, 0)]
                    univ_ids = words
                else:
                    pairs = re.findall(r'-?\d+\s*:\s*-?\d+', fill)
                    i_colon = fill.rfind(':')
                    univ_ids = fill[i_colon + 1:].split()[1:]

                    if not pairs:
                        raise ValueError('Cant find lattice specification')

                    ranges = [tuple(map(int, pairs[i].split(':')))
                              for i in range(3)]
                    for lower, upper in ranges:
                        assert upper >= lower
                univ_ids = np.asarray(univ_ids, dtype=int)

                # A finite lattice with a single axial layer becomes a 2D
                # lattice whose universes are translated to the layer, so that
                # the top and bottom of the layer are not lattice surfaces
                # coincident with those of the cell that contains it
                k1, k2 = ranges[2]
                if len(vectors) == (4 if hexagonal else 3) and k1 == k2 \
                        and not inf_lattice:
                    center[2] = -k1*vectors[-1][2]
                    vectors = vectors[:-1]

                # Check for universe ID same as the ID assigned to the cell
                # itself -- since OpenMC can't handle this directly, we need
                # to create an extra cell/universe to fill in the lattice. The
                # same universe fills the hexagons that are not part of the
                # FILL array of a hexagonal lattice.
                if hexagonal or np.any(univ_ids == uid):
                    extra_cell = openmc.Cell(
                        fill=mat if cell_material_id > 0 else None)
                    filler = openmc.Universe(cells=[extra_cell])
                    univ_ids[univ_ids == uid] = filler.id

                    # Put it in universes dictionary so that get_universe
                    # works correctly
                    universes[filler.id] = filler

                # If center of MCNP lattice element is not (0,0,0), we need
                # to translate the universe
                if not np.all(center == 0.0):
                    for uid in np.unique(univ_ids):
                        # Create translated universe
                        trans_cell = openmc.Cell(fill=get_universe(uid))
                        trans_cell.translation = -center
                        u = openmc.Universe(cells=[trans_cell])
                        universes[u.id] = u

                        # Replace original universes with translated ones
                        univ_ids[univ_ids == uid] = u.id

                if hexagonal:
                    _fill_hex_lattice(lattice, center, vectors, ranges,
                                      univ_ids, get_universe, filler)
                else:
                    _fill_rect_lattice(lattice, center, vectors, ranges,
                                       univ_ids, get_universe)

                # For infinite lattices, set the outer universe
                if inf_lattice:
                    lattice.outer = get_universe(univ_ids[0])

                cell._lattice = True
            else:
                # Cell filled with universes
                if 'fill' in c['parameters']:
                    uid, ftrans = _CELL_FILL_RE.search(c['parameters']['fill']).groups()
                    use_degrees = False
                else:
                    uid, ftrans = _CELL_FILL_RE.search(c['parameters']['*fill']).groups()
                    use_degrees = True

                # First assign fill based on whether it is a universe/lattice
                uid = int(uid)
                if uid not in universes:
                    for ci in cells:
                        if 'u' in ci['parameters']:
                            if abs(int(ci['parameters']['u'])) == uid:
                                if 'lat' in ci['parameters']:
                                    universes[uid] = _new_lattice(
                                        ci['parameters'], uid)
                                else:
                                    universes[uid] = openmc.Universe(uid)
                                break
                cell.fill = universes[uid]

                # Set fill transformation
                if ftrans is not None:
                    ftrans = ftrans.split()
                    if len(ftrans) > 3:
                        vector = np.array([float(x) for x in ftrans[:3]])
                        if len(ftrans) == 13:
                            if int(ftrans[12]) == -1:
                                vector *= -1

                        cell.translation = tuple(vector)
                        rotation_matrix = np.array([float(x) for x in ftrans[3:12]]).reshape((3, 3))
                        if use_degrees:
                            rotation_matrix = np.cos(rotation_matrix * pi/180.0)
                        cell.rotation = rotation_matrix
                    elif len(ftrans) < 3:
                        assert len(ftrans) == 1
                        tr_num = int(ftrans[0])
                        translation, rotation = data['tr'][tr_num]
                        cell.translation = translation
                        if rotation is not None:
                            cell.rotation = rotation.T
                    else:
                        cell.translation = tuple(float(x) for x in ftrans)

        elif c['material'] > 0:
            cell.fill = mat

        if 'vol' in c["parameters"]:
            cell.volume = float(c["parameters"]["vol"])

        if not hasattr(cell, '_lattice'):
            openmc_cells[c['id']] = cell

    # Expand shorthand notation
    def replace_complement(region, cells):
        if isinstance(region, (openmc.Intersection, openmc.Union)):
            for n in region:
                replace_complement(n, cells)
        elif isinstance(region, openmc.Complement):
            if isinstance(region.node, openmc.Halfspace):
                region.node = cells[region.node.surface.id].region

    for cell in openmc_cells.values():
        replace_complement(cell.region, openmc_cells)
    return universes


def mcnp_to_model(filename, merge_surfaces: bool = True, expand_elements: bool = True) -> openmc.Model:
    """Convert MCNP input to OpenMC model

    Parameters
    ----------
    filename : str
        Path to MCNP file
    merge_surfaces : bool
        Whether to remove redundant surfaces when the geometry is exported.

    Returns
    -------
    openmc.Model
        Equivalent OpenMC model

    """

    cells, surfaces, data = parse(filename)

    openmc_materials = get_openmc_materials(data['materials'], expand_elements)
    openmc_surfaces = get_openmc_surfaces(surfaces, data)
    openmc_universes = get_openmc_universes(cells, openmc_surfaces,
                                            openmc_materials, data)

    geometry = openmc.Geometry(openmc_universes[0])
    geometry.merge_surfaces = merge_surfaces
    materials = openmc.Materials(geometry.get_all_materials().values())

    settings = openmc.Settings()
    settings.batches = 40
    settings.inactive = 20
    settings.particles = 100
    settings.output = {'summary': True}

    # Determine bounding box for geometry
    all_volume = openmc.Union([cell.region for cell in
                                geometry.root_universe.cells.values()])
    ll, ur = all_volume.bounding_box
    src_class = getattr(openmc, 'IndependentSource')
    if src_class is None:
        src_class = openmc.Source
    if np.any(np.isinf(ll)) or np.any(np.isinf(ur)):
        settings.source = src_class(space=openmc.stats.Point())
    else:
        settings.source = src_class(space=openmc.stats.Point((ll + ur)/2))

    return openmc.Model(geometry, materials, settings)


def mcnp_str_to_model(text: str, **kwargs):
    # Write string to a temporary file
    with tempfile.NamedTemporaryFile('w', delete=False) as fp:
        fp.write(text)

    # Parse model from file
    model = mcnp_to_model(fp.name, **kwargs)

    # Remove temporary file and return model
    os.remove(fp.name)
    return model


def mcnp_to_openmc():
    """Command-line interface for converting MCNP model"""
    parser = argparse.ArgumentParser()
    parser.add_argument('mcnp_filename')
    parser.add_argument('--merge-surfaces', action='store_true',
                        help='Remove redundant surfaces when exporting XML')
    parser.add_argument('--no-merge-surfaces', dest='merge_surfaces', action='store_false',
                        help='Do not remove redundant surfaces when exporting XML')
    parser.add_argument('--expand-elements', action='store_true',
                        help='Expand elements to their constituent isotopes')
    parser.add_argument('--no-expand-elements', dest='expand_elements', action='store_false',
                        help='Do not expand elements to their constituent isotopes')
    parser.add_argument('-o', '--output', default='model.xml',
                        help='Name for the OpenMC model XML file')
    parser.add_argument('-s', '--separate-xml', action='store_true',
                        help='Write separate XML files')
    parser.set_defaults(merge_surfaces=True)
    parser.set_defaults(expand_elements=True)
    args = parser.parse_args()

    model = mcnp_to_model(args.mcnp_filename, args.merge_surfaces, args.expand_elements)
    if args.separate_xml:
        model.export_to_xml()
    else:
        model.export_to_model_xml(args.output)
