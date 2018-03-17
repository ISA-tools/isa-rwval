#!/usr/bin/env python

import re
_RX_PARAMETER_VALUE = re.compile('Parameter Value\[(.*?)\]')
_RX_FACTOR_VALUE = re.compile('Factor Value\[(.*?)\]')
_RX_SOURCE = re.compile('^Source Name$')
_RX_SAMPLE = re.compile('^Sample Name$')
_RX_CHARACTERISTICS = re.compile('^Characteristics\[(.*?)\]$')
_RX_COMMENT = re.compile('^Comment\[(.*?)\]$')
_RX_UNIT = re.compile('^Unit$')
_RX_TERM_SOURCE_REF = re.compile('^Term Source REF$')
_RX_TERM_ACCESSION_NUMBER = re.compile('^Term Accession Number$')
_RX_PROTOCOL_REF = re.compile('^Protocol REF$')

import csv
import itertools
import networkx as nx

def _is_node(x):
    is_node = False
    if _RX_SOURCE.match(x):
        is_node = True
    elif _RX_SAMPLE.match(x):
        is_node = True
    elif _RX_PROTOCOL_REF.match(x):
        is_node = True
    return is_node


def _pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def _parse_material(header, row):
    from isatools.model import Source, Sample, Extract, LabeledExtract
    obj = None
    node_id = row[0]
    obj_label = header[0]
    characteristics = []
    factor_values = []
    comments = []
    if obj_label == 'Source Name':
        obj = Source(name=node_id)
        obj.characteristics = characteristics
        obj.comments = comments
    elif obj_label == 'Sample Name':
        obj = Sample(name=node_id)
        obj.characteristics = characteristics
        obj.comments = comments
        obj.factor_values = factor_values
    elif obj_label == 'Sample Name':
        obj = Sample(name=node_id)
        obj.characteristics = characteristics
        obj.comments = comments
        obj.factor_values = factor_values
    elif obj_label == 'Extract Name':
        obj = Extract(name=node_id)
        obj.characteristics = characteristics
        obj.comments = comments
    elif obj_label == 'Labeled Extract Name':
        obj = LabeledExtract(name=node_id)
        obj.characteristics = characteristics
        obj.comments = comments
    return obj


def _parse_process(header, row):
    from isatools.model import Process, Protocol
    obj = None
    parameter_values = []
    comments = []
    obj = Process()
    if header[0] == 'Protocol REF':
        protocol_key = row[0]
        obj.executes_protocol = Protocol(name=protocol_key)
    obj.parameter_values = parameter_values
    obj.comments = comments
    if 'Assay Name' in header:
        process_name = row[header.index('Assay Name')]
        obj.name = process_name
    return obj


with open('/Users/dj/Development/ISA/isatools-core/tests/data/tab/BII-I-1/s_BII-S-1.txt') as fp:
    reader = csv.reader(fp, delimiter='\t')
    header = next(reader)  # always assume the first row is the header
    print('Header is: {header}'.format(header=header))
    len_header = len(header)
    print('Length of header is {len_header}'.format(len_header=len(header)))
    G_list = []
    for rn, row in enumerate(reader):
        if len(row) != len_header:
            print('Skipping row {rn} as row length does not match header length'.format(rn=rn))  # warn if a row is skipped
        else:
            G = nx.DiGraph()
            indices = [i for i, x in enumerate(header) if _is_node(x)]
            if not _is_node(header[-1]):
                indices.append(-1)
            L, R = None, None
            for x, y in _pairwise(indices):
                print(header[x:y])
                if header[x] in ('Source Name', 'Sample Name'):
                    obj = _parse_material(header[x:y], row[x:y])
                    print('popped ', obj)
                    print('L is {} and R is {}'.format(L, R))
                    if not L and not R:
                        L = obj
                    elif L and not R:
                        R = obj
                    else:
                        G.add_edge(L, R)
                        L, R = R, None
                elif header[x] == 'Protocol REF':
                    obj = _parse_process(header[x:y], row[x:y])
                    print('popped ', obj)
                    print('L is {} and R is {}'.format(L, R))
                    if not L and not R:
                        L = obj
                    elif L and not R:
                        R = obj
                    else:
                        G.add_edge(L, R)
                        L, R = R, None
            G_list.append(G.edges)
    print(G_list)
