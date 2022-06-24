# SPDX-FileCopyrightText: 2022 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

from collections import defaultdict
from copy import deepcopy
from math import pi
import re

import numpy as np


_CELL1_RE = re.compile(r'\s*(\d+)\s+(\d+)([ \t0-9:#().dDeE\+-]+)((?:\S+\s*=.*)?)')
_CELL_KEYWORDS_RE = re.compile(r'[A-Za-z: \t]+=[-0-9:() \t]+')
_CELL_FILL_RE = re.compile(r'\s*(\d+)\s*(?:\((.*)\))?')
_CELL2_RE = re.compile(r'\s*(\d+)\s+like\s+(\d+)\s+but\s+(\S+\s*=.*)')
_SURFACE_RE = re.compile(r'\s*(\*?\d+)(\s*[-0-9]+)?\s+(\S+)((?:\s+\S+)+)')
_MATERIAL_RE = re.compile(r'\s*[Mm](\d+)((?:\s+\S+)+)')
_TR_RE = re.compile(r'\s*(\*)?[Tt][Rr](\d+)\s+(.*)')
_SAB_RE = re.compile(r'\s*[Mm][Tt](\d+)((?:\s+\S+)+)')
_MODE_RE = re.compile(r'\s*mode(?:\s+\S+)*')
_COMPLEMENT_RE = re.compile(r'(#)(\d+)')
_REPEAT_RE = re.compile(r'(\d+)\s+(\d+)[rR]')
_NUM_RE = re.compile(r'(\d)([+-])(\d)')


def float_(val):
    """Convert scientific notation literals that don't have an 'e' in them to float"""
    return float(_NUM_RE.sub(r'\1e\2\3', val))


def parse_cell(line):
    """Parse cell card into dictionary of information

    Parameters
    ----------
    line : str
        Single MCNP cell card

    Returns
    -------
    dict
        Dictionary with cell information

    """
    # Helper function for population parameters dictionary
    def get_parameters(s, params):
        # TODO: Use _CELL_KEYWORDS_RE?
        words = (s + ' after').split('=')
        for before, after in zip(words[:-1], words[1:]):
            key = before.split()[-1]
            value = ' '.join(after.split()[:-1])
            params[key] = value

    if 'like' in line.lower():
        # Handle LIKE n BUT form
        m = _CELL2_RE.match(line.lower())
        if m is not None:
            g = m.groups()
            parameters = {}
            if len(g) == 3:
                get_parameters(g[-1], parameters)
            return {'id': int(g[0]), 'like': int(g[1]), 'parameters': parameters}

    else:
        # Handle normal form
        m = _CELL1_RE.match(line.lower())
        if m is not None:
            g = m.groups()
            if g[1] == '0':
                density = None
                region = g[2].strip()
            else:
                words = g[2].split()
                density = float_(words[0])
                region = ' '.join(words[1:])
            parameters = {}
            if g[-1]:
                get_parameters(g[-1], parameters)
            return {'id': int(g[0]), 'material': int(g[1]), 'density': density,
                    'region': region, 'parameters': parameters}


def resolve_likenbut(cells):
    """Resolve LIKE n BUT by copying attributes from like cell

    Parameters
    ----------
    cells : list of dict
        Cell information for each cell

    """
    # Create dictionary mapping ID to cell dict
    cell_by_id = {c['id']: c for c in cells}
    for cell in cells:
        if 'like' in cell:
            cell_id = cell['id']
            like_cell_id = cell['like']
            params = cell['parameters']

            # Clear dictionary and copy information from like cell
            cell.clear()
            cell.update(deepcopy(cell_by_id[like_cell_id]))

            # Update ID and specified parameters
            cell['id'] = cell_id
            for key, value in params.items():
                if key == 'mat':
                    cell['material'] = int(value)
                elif key == 'rho':
                    cell['density'] = float_(value)
                else:
                    cell['parameters'][key] = value


def parse_surface(line):
    """Parse surface card into dictionary of information

    Parameters
    ----------
    line : str
        Single MCNP surface card

    Returns
    -------
    dict
        Dictionary with surface information

    """
    m = _SURFACE_RE.match(line)
    if m is not None:
        g = m.groups()
        surface = {}
        if '*' in g[0]:
            surface['reflective'] = True
            uid = int(g[0][1:])
        else:
            surface['reflective'] = False
            uid = int(g[0])
        surface.update({'id': uid, 'mnemonic': g[2].lower(),
                        'coefficients': g[3].strip()})
        if g[1] is not None:
            if int(g[1]) < 0:
                surface['periodic'] = int(g[1])
                # TODO: Move into OpenMC conversion
                raise NotImplementedError('Periodic boundary conditions not supported')
            else:
                surface['tr'] = int(g[1])
        return surface
    else:
        raise NotImplementedError("Unable to convert surface card: {}".format(line))


def parse_data(section):
    """Parse data block into dictionary of information

    Parameters
    ----------
    line : str
        MCNP data block

    Returns
    -------
    dict
        Dictionary with data-block information

    """

    data = {'materials': defaultdict(dict), 'tr': {}}

    lines = section.split('\n')
    for line in lines:
        if _MATERIAL_RE.match(line):
            g = _MATERIAL_RE.match(line).groups()
            spec = g[1].split()
            try:
                nuclides = list(zip(spec[::2], map(float_, spec[1::2])))
            except Exception:
                raise ValueError('Invalid material specification?')
            uid = int(g[0])
            data['materials'][uid].update({'id': uid, 'nuclides': nuclides})
        elif _SAB_RE.match(line):
            g = _SAB_RE.match(line).groups()
            uid = int(g[0])
            spec = g[1].split()
            data['materials'][uid]['sab'] = spec
        elif line.lower().startswith('mode'):
            words = line.split()
            data['mode'] = words[1:]
        elif line.lower().startswith('kcode'):
            words = line.split()
            data['kcode'] = {'n_particles': int(words[1]),
                             'initial_k': float_(words[2]),
                             'inactive': int(words[3]),
                             'batches': int(words[4])}
        elif _TR_RE.match(line):
            g = _TR_RE.match(line).groups()
            use_degrees = g[0] is not None
            tr_num = int(g[1])
            values = g[2].split()
            if len(values) >= 3:
                displacement = np.array([float(x) for x in values[:3]])
            if len(values) >= 12:
                rotation = np.array([float(x) for x in values[3:12]]).reshape((3,3)).T
                if use_degrees:
                    rotation = np.cos(rotation * pi/180.0)
            else:
                rotation = None
            data['tr'][tr_num] = (displacement, rotation)

    return data


def split_mcnp(filename):
    """Split MCNP file into three strings, one for each block

    Parameters
    ----------
    filename : str
        Path to MCNP file

    Returns
    -------
    list of str
        List containing one string for each block

    """
    # Find beginning of cell section
    text = open(filename, 'r').read()
    m = re.search(r'^[ \t]*(\d+)[ \t]+', text, flags=re.MULTILINE)
    text = text[m.start():]
    return re.split('\n[ \t]*\n', text)


def sanitize(section):
    """Sanitize one section of an MCNP input

    This function will remove comments, join continuation lines into a single
    line, and expand repeated numbers explicitly.

    Parameters
    ----------
    section : str
        String representing one section of an MCNP input

    Returns
    -------
    str
        Sanitized input section

    """

    # Remove end-of-line comments
    section = re.sub('\$.*$', '', section, flags=re.MULTILINE)

    # Remove comment cards
    section = re.sub('^[ \t]*?[cC].*?$\n?', '', section, flags=re.MULTILINE)

    # Turn continuation lines into single line
    section = re.sub('&.*\n', ' ', section)
    section = re.sub('\n {5}', ' ', section)

    # Expand repeated numbers
    m = _REPEAT_RE.search(section)
    while m is not None:
        section = _REPEAT_RE.sub(' '.join((int(m.group(2)) + 1)*[m.group(1)]),
                                 section, 1)
        m = _REPEAT_RE.search(section)
    return section


def parse(filename):
    """Parse an MCNP file and return information from three main blocks

    Parameters
    ----------
    filename : str
        Path to MCNP file

    Returns
    -------
    cells : list
        List of dictionaries, where each dictionary contains information for one cell
    surfaces : list
        List of dictionaries, where each dictionary contains information for one surface
    data : dict
        Dictionary containing data-block information, including materials

    """
    # Split file into main three sections (cells, surfaces, data)
    sections = split_mcnp(filename)

    # Sanitize lines (remove comments, continuation lines, etc.)
    cell_section = sanitize(sections[0])
    surface_section = sanitize(sections[1])
    data_section = sanitize(sections[2])

    cells = [parse_cell(x) for x in cell_section.strip().split('\n')]
    surfaces = [parse_surface(x) for x in surface_section.strip().split('\n')]
    data = parse_data(data_section)

    # Replace LIKE n BUT with actual parameters
    resolve_likenbut(cells)

    return cells, surfaces, data
