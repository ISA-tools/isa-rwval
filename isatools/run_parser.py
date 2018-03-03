#!/usr/bin/env python

import csv
import networkx as nx

import re

_RX_DATA = re.compile('data\[(.*?)\]')
_RX_DOI = re.compile('(10[.][0-9]{4,}(?:[.][0-9]+)*/(?:(?![%"#? ])\\S)+)')
_RX_PMID = re.compile('[0-9]{8}')
_RX_PMCID = re.compile('PMC[0-9]{8}')
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

TOKEN_MAP = {
    'Source Name': _RX_SOURCE,
    'Sample Name': _RX_SAMPLE,
    'Characteristics':_RX_CHARACTERISTICS,
    'Comment': _RX_COMMENT,
    'Protocol REF': _RX_PROTOCOL_REF
}

from model import Sample
from model import Source

CLASS_MAP = {
    'Source Name': Source,
    'Sample Name': Sample
}


def checkfor(token, nextcell):
    if isinstance(token, str):
        token_regex = TOKEN_MAP[token]
        if not token_regex.match(nextcell[0]):
            raise SyntaxError(
                'Next token must begin with {expected} but got {actual}'
                .format(expected=token, actual=nextcell[0]))


def parse_characteristics(nextcell, mr_iter, cn):
    from model import Characteristic, OntologyAnnotation, OntologySource
    category_term = _RX_CHARACTERISTICS.findall(nextcell[0])[-1]
    category = OntologyAnnotation(term=category_term)
    characteristic = Characteristic(category=category)
    value = nextcell[1]
    characteristic.value = value
    nextcell = next(mr_iter)
    cn += 1
    if _RX_TERM_SOURCE_REF.match(nextcell[0]):
        value = OntologyAnnotation(term=value)
        value.term_source = OntologySource(name=nextcell[1])
        nextcell = next(mr_iter)
        cn += 1
        checkfor(_RX_TERM_ACCESSION_NUMBER, nextcell)
        value.term_accession = nextcell[1]
        characteristic.value = value
        nextcell = next(mr_iter)
        cn += 1
    return nextcell, mr_iter, characteristic, cn


def parse_factor_values(nextcell, mr_iter, cn):
    from model import FactorValue, StudyFactor, OntologyAnnotation, OntologySource
    factor_name = _RX_FACTOR_VALUE.findall(nextcell[0])[-1]
    factor = StudyFactor(name=factor_name)
    factor_value = FactorValue(factor_name=factor_name)
    value = nextcell[1]
    factor_value.value = value
    nextcell = next(mr_iter)
    cn += 1
    if _RX_TERM_SOURCE_REF.match(nextcell[0]):
        value = OntologyAnnotation(term=value)
        value.term_source = OntologySource(name=nextcell[1])
        nextcell = next(mr_iter)
        cn += 1
        checkfor(_RX_TERM_ACCESSION_NUMBER, nextcell)
        value.term_accession = nextcell[1]
        factor_value.value = value
        try:
            nextcell = next(mr_iter)
            cn += 1
        except StopIteration:
            nextcell = None # EOF
    elif _RX_UNIT.match(nextcell[0]):
        # checkfor(_RX_UNIT, nextcell)
        unit = OntologyAnnotation(term=nextcell[1])
        nextcell = next(mr_iter)
        cn += 1
        #checkfor(_RX_TERM_SOURCE_REF, nextcell)
        unit.term_source = OntologySource(name=nextcell[1])
        nextcell = next(mr_iter)
        cn += 1
        checkfor(_RX_TERM_ACCESSION_NUMBER, nextcell)
        unit.term_accession = nextcell[1]
        factor_value.unit = unit
        try:
            nextcell = next(mr_iter)
            cn += 1
        except StopIteration:
            nextcell = None # EOF
    return nextcell, mr_iter, factor_value, cn


def parse_comment(nextcell, mr_iter, cn):
    from model import Comment
    name = _RX_COMMENT.findall(nextcell[0])[-1]
    comment = Comment(name=name)
    value = nextcell[1]
    comment.value = value
    nextcell = next(mr_iter)
    cn += 1
    return nextcell, mr_iter, comment, cn


def parse_material(context, nextcell, mr_iter, rn, cn):
    if not nextcell[1].strip():
        raise SyntaxError("'{context}' value must not be empty".format(
            context=context))
    else:
        materialClass = CLASS_MAP[context]
        print('parsing {context} at {loc}'.format(context=context, loc=(rn, cn)))
        material = materialClass(name=nextcell[1])
        nextcell = next(mr_iter)
        cn += 1
        while _RX_CHARACTERISTICS.match(nextcell[0]) or _RX_FACTOR_VALUE.match(nextcell[0]) or _RX_COMMENT.match(nextcell[0]):
            if _RX_CHARACTERISTICS.match(nextcell[0]):
                print('parsing characteristic at {loc}'.format(loc=(rn, cn)))
                nextcell, mr_iter, characteristic, cn = parse_characteristics(nextcell, mr_iter, cn)
                material.characteristics.append(characteristic)
            elif _RX_FACTOR_VALUE.match(nextcell[0]):
                print('parsing factor value at {loc}'.format(loc=(rn, cn)))
                nextcell, mr_iter, factor_value, cn = parse_factor_values(nextcell, mr_iter, cn)
                material.factor_values.append(factor_value)
            elif _RX_COMMENT.match(nextcell[0]):
                print('parsing comment at {loc}'.format(loc=(rn, cn)))
                nextcell, mr_iter, comment, cn = parse_comment(nextcell, mr_iter, cn)
                material.comments.append(comment)
    return nextcell, mr_iter, material, cn


def parse_protocol_ref(context, nextcell, mr_iter, rn, cn):
    if not nextcell[1].strip():
        raise SyntaxError("'{context}' value must not be empty".format(
            context=context))
    else:
        from model import Process, Protocol
        print('parsing process (Protocol REF) at {loc}'.format(loc=(rn, cn)))
        process = Process(executes_protocol=Protocol(name=nextcell[1]))
        nextcell = next(mr_iter)
        cn += 1
    return nextcell, mr_iter, process, cn


def parse_row(header, row, rn):
    #keygen = lambda x: ('.'.join(x) if x[1] else None)
    G = nx.DiGraph()
    mr_iter = iter(map(lambda x: (x[0], x[1]), zip(header, row)))
    nextcell = next(mr_iter)
    cn = 0
    checkfor('Source Name', nextcell)
    nextcell, mr_iter, M, cn = parse_material('Source Name', nextcell, mr_iter, rn, cn)
    while nextcell:
        checkfor('Protocol REF', nextcell)
        print('Protocol REF', nextcell, mr_iter, rn, cn)
        nextcell, mr_iter, P, cn = parse_protocol_ref('Protocol REF', nextcell, mr_iter, rn, cn)
        G.add_edge(M, P)
        checkfor('Sample Name', nextcell)
        nextcell, mr_iter, M, cn = parse_material('Sample Name', nextcell, mr_iter, rn, cn)
        G.add_edge(P, M)
    return G


with open('/Users/dj/Development/ISA/isatools-core/tests/data/tab/BII-I-1/s_BII-S-1.txt') as fp:
    reader = csv.reader(fp, delimiter='\t')
    header = next(reader)  # always assume the first row is the header
    print('Header is: {header}'.format(header=header))
    len_header = len(header)
    print('Length of header is {len_header}'.format(len_header=len(header)))
    list_of_digraphs = []
    for rn, row in enumerate(reader):
        G = parse_row(header, row, rn)
        print(G.nodes)

"""
with open('/Users/dj/Development/ISA/isatools-core/tests/data/tab/BII-I-1/s_BII-S-1.txt') as fp:
    reader = csv.reader(fp, delimiter='\t')
    header = next(reader)  # always assume the first row is the header
    print('Header is: {header}'.format(header=header))
    len_header = len(header)
    print('Length of header is {len_header}'.format(len_header=len(header)))
    list_of_digraphs = []
    for rn, row in enumerate(reader):
        if len(row) != len_header:
            print('Skipping row {rn} as row length does not match header length'.format(rn=rn)) # warn if a row is skipped
        else:
            print('Row {rn} is: {row_content}'.format(rn=rn, row_content=row))
            mapped_row = map(
                lambda x: ('.'.join(x) if x[1] else None, x[0], x[1]), zip(header, row))
            print('Header-Row KVs for row {rn}: {row_as_kvs}'.format(
                rn=rn, row_as_kvs=mapped_row))
            g = nx.DiGraph()
            filtered_row = []
            for item in mapped_row:
                if item[0] and 'Name' in item[1] or 'Protocol' in item[1]:
                    filtered_row.append(item)
            print('Filtered KVs for row {rn}: {filtered_row}'.format(
                rn=rn, filtered_row=filtered_row))
            from itertools import tee, izip

            def pairwise(iterable):
                a, b = tee(iterable)
                next(b, None)
                return izip(a, b)

            for u, v in pairwise(filtered_row):
                cn1 = header.index(u[1])
                cn2 = header.index(v[1])
                print('Adding edge for row {rn} and objects found at column '
                      '{cn1} and {cn2}'.format(rn=rn, cn1=cn1, cn2=cn2))
                g.add_edge(u[0], v[0])
            print(g.edges)
"""

"""
from .model import *
class StudyGraph(nx.DiGraph):

    def add_isa_edge(self, u, v):

        if isinstance(u, (Source, Sample)):
            assert isinstance(v, Protocol)
        if isinstance(u, Protocol):
            assert isinstance(v, Sample)
        self.add_edge(u, v)
"""