from collections import defaultdict
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
    return float(_NUM_RE.sub(r'\1e\2\3', val))


def parse_cell(line):
    if 'like' in line.lower():
        raise NotImplementedError('like N but form not supported')
        # Handle LIKE n BUT form
        m = _CELL2_RE.match(line.lower())
        if m is not None:
            g = m.groups()
            parameters = (dict(kv.split('=') for kv in g[-1].split())
                          if len(g) == 3 else {})
            return {'id': g[0], 'like': g[1], 'parameters': parameters}

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
                words = (g[-1] + ' after').split('=')
                for before, after in zip(words[:-1], words[1:]):
                    key = before.split()[-1]
                    value = ' '.join(after.split()[:-1])
                    parameters[key] = value
            return {'id': int(g[0]), 'material': int(g[1]), 'density': density,
                    'region': region, 'parameters': parameters}


def parse_surface(line):
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
                raise NotImplementedError('Periodic boundary conditions not supported')
            else:
                surface['tr'] = int(g[1])
        return surface
    else:
        raise NotImplementedError("Unable to convert surface card: {}".format(line))


def parse_data(section):
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


def parse(filename):
    # Find beginning of cell section
    text = open(filename, 'r').read()
    m = re.search(r'^[ \t]*(\d+)[ \t]+', text, flags=re.MULTILINE)
    text = text[m.start():]

    sections = re.split('\n[ \t]*\n', text)
    for i in range(len(sections)):
        # Remove end-of-line comments
        sections[i] = re.sub('\$.*$', '', sections[i], flags=re.MULTILINE)

        # Remove comment cards
        sections[i] = re.sub('^[ \t]*?[cC].*?$\n?', '', sections[i], flags=re.MULTILINE)

        # Turn continuation lines into single line
        sections[i] = re.sub('&.*\n', ' ', sections[i])
        sections[i] = re.sub('\n {5}', ' ', sections[i])

        # Expand repeated numbers
        m = _REPEAT_RE.search(sections[i])
        while m is not None:
            sections[i] = _REPEAT_RE.sub(' '.join((int(m.group(2)) + 1)*[m.group(1)]),
                                         sections[i], 1)
            m = _REPEAT_RE.search(sections[i])

    cells = list(map(parse_cell, sections[0].strip().split('\n')))
    surfaces = list(map(parse_surface, sections[1].strip().split('\n')))
    data = parse_data(sections[2])

    return cells, surfaces, data


