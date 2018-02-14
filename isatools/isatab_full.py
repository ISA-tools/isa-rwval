
from __future__ import absolute_import

from bisect import bisect_left, bisect_right
from collections import OrderedDict

import logging
import os
import re

import numpy as np
import pandas as pd

from isatools.isatab_meta import AbstractParser, InvestigationParser
from isatools.model import *

log = logging.getLogger(__name__)

class TableParser(AbstractParser):

    DATA_FILE_LABELS = (
        'Raw Data File', 'Derived Spectral Data File',
        'Derived Array Data File', 'Array Data File',
        'Protein Assignment File', 'Peptide Assignment File',
        'Post Translational Modification Assignment File',
        'Acquisition Parameter Data File', 'Free Induction Decay Data File',
        'Derived Array Data Matrix File', 'Image File', 'Derived Data File',
        'Metabolite Assignment File', 'Raw Spectral Data File')
    MATERIAL_LABELS = ('Source Name', 'Sample Name', 'Extract Name',
                       'Labeled Extract Name')
    OTHER_MATERIAL_LABELS = ('Extract Name', 'Labeled Extract Name')
    NODE_LABELS = DATA_FILE_LABELS + MATERIAL_LABELS + OTHER_MATERIAL_LABELS
    ASSAY_LABELS = ('Assay Name', 'MS Assay Name', 'Hybridization Assay Name',
                    'Scan Name', 'Data Transformation Name',
                    'Normalization Name')
    ALL_LABELS = NODE_LABELS + ASSAY_LABELS + tuple(['Protocol REF'])

    # REGEXES
    RX_I_FILE_NAME = re.compile(r'i_(.*?)\.txt')
    RX_DATA = re.compile(r'data\[(.*?)\]')
    RX_COMMENT = re.compile(r'Comment\[(.*?)\]')
    RX_DOI = re.compile(r'(10[.][0-9]{4,}(?:[.][0-9]+)*/(?:(?![%"#? ])\\S)+)')
    RX_PMID = re.compile(r'[0-9]{8}')
    RX_PMCID = re.compile(r'PMC[0-9]{8}')
    RX_CHARACTERISTICS = re.compile(r'Characteristics\[(.*?)\]')
    RX_PARAMETER_VALUE = re.compile(r'Parameter Value\[(.*?)\]')
    RX_FACTOR_VALUE = re.compile(r'Factor Value\[(.*?)\]')
    RX_INDEXED_COL = re.compile(r'(.*?)\.\d+')

    @staticmethod
    def _find_lt(a, x):
        i = bisect_left(a, x)
        if i:
            return a[i - 1]
        else:
            return -1

    @staticmethod
    def _find_gt(a, x):
        i = bisect_right(a, x)
        if i != len(a):
            return a[i]
        else:
            return -1

    @staticmethod
    def _clean_label(label):
        for clean_label in TableParser.ALL_LABELS:
            if label.startswith(clean_label):
                return clean_label

    def __init__(self):
        self.node_map = OrderedDict()
        self.process_map = OrderedDict()
        self.ontology_sources = OrderedDict()
        self.characteristic_categories = OrderedDict()
        self.unit_categories = OrderedDict()

    def _insert_missing_protocol_refs(self, df):
        raise NotImplementedError()  # TODO: Finish implementation
        df = df[[x for x in df.columns if
                 x.startswith(TableParser.ALL_LABELS)]]
        labels = df.columns
        nodes_index = [i for i, x in enumerate(labels) if
                       x in TableParser.NODE_LABELS]
        missing_protocol_ref_indicies = []
        for cindex, label in enumerate(labels):
            if not label.startswith('Protocol REF'):
                output_node_index = self._find_gt(nodes_index, cindex)
                if output_node_index > -1:
                    output_node_label = labels[output_node_index]
                    if output_node_label in TableParser.MATERIAL_LABELS + \
                        TableParser.DATA_FILE_LABELS + TableParser.ASSAY_LABELS:
                        missing_protocol_ref_indicies.append(output_node_index)
        offset = 0
        for i in reversed(missing_protocol_ref_indicies):
            inferred_protocol_type = ''
            leftcol = labels[self._find_lt(nodes_index, i)]
            rightcol = labels[i]
            if leftcol == 'Source Name' and rightcol == 'Sample Name':
                inferred_protocol_type = 'sample collection'
            elif leftcol == 'Sample Name' and rightcol == 'Extract Name':
                inferred_protocol_type = 'extraction'
            elif leftcol == 'Extract Name' and rightcol == \
                    'Labeled Extract Name':
                inferred_protocol_type = 'labeling'
            elif leftcol == 'Labeled Extract Name' and rightcol in (
                    'Assay Name', 'MS Assay Name'):
                inferred_protocol_type = 'library sequencing'
            elif leftcol == 'Extract Name' and rightcol in (
                    'Assay Name', 'MS Assay Name'):
                inferred_protocol_type = 'library preparation'
            elif leftcol == 'Scan Name' and rightcol == 'Raw Data File':
                inferred_protocol_type = 'data acquisition'
            elif leftcol == 'Assay Name' and rightcol == 'Normalization Name':
                inferred_protocol_type = 'normalization'
            elif leftcol == 'Normalization Name' and \
                            rightcol == 'Data Transformation Name':
                inferred_protocol_type = 'data transformation'
            elif leftcol == 'Raw Data File' and \
                            rightcol == 'Metabolite Identification File':
                inferred_protocol_type = 'metabolite identification'
            elif leftcol == 'Raw Data File' and \
                            rightcol == 'Protein Identification File':
                inferred_protocol_type = 'metabolite identification'
            # Force use of unknown protocol always, until we can insert missing
            # protocol from above inferences into study metadata
            log.info('Inserting protocol %s in between %s and %s', 
                inferred_protocol_type if inferred_protocol_type != '' else 'unknown',
                leftcol, rightcol)
            protocol_ref_cols = [x for x in labels if
                                 x.startswith('Protocol REF')]
            num_protocol_refs = len(protocol_ref_cols)
            df.insert(i, 'Protocol REF.{}'.format(num_protocol_refs + offset),
                      'unknown' if inferred_protocol_type == '' else
                      inferred_protocol_type)
            offset += 1
        return df

    def _attach_factors(self, df):
        df = df [[x for x in df.columns if x.startswith(
            'Sample Name', 'Factor Value')]]
        for _, row in df.iterrows():
            sample = self.node_map['Sample Name.{sample_name}'.format(
                sample_name=row['Sample Name'])]
            log.info('attaching factor to sample: %s', sample.name)

    def _make_process_sequence(self, df):
        """This function builds the process sequences and links nodes to
        processes based on node keys calculated in the ISA-Tab tables"""
        df = df[[x for x in df.columns if
                 x.startswith(TableParser.ALL_LABELS)]]
        process_key_sequences = []
        for _, row in df.iterrows():
            process_key_sequence = []
            labels = df.columns
            nodes_index = [i for i, x in enumerate(labels) if
                           x in TableParser.NODE_LABELS]
            for cindex, label in enumerate(labels):
                val = row[label]
                if label.startswith('Protocol REF') and val != '':
                    output_node_index = self._find_gt(nodes_index, cindex)
                    if output_node_index > -1:
                        output_node_label = labels[output_node_index]
                        output_node_val = row[output_node_label]
                    input_node_index = self._find_lt(nodes_index, cindex)
                    if input_node_index > -1:
                        input_node_label = labels[input_node_index]
                        input_node_val = row[input_node_label]
                    input_nodes_with_prot_keys = df.loc[
                        df[labels[cindex]] == val].groupby(
                        [labels[cindex], labels[input_node_index]]).size()
                    output_nodes_with_prot_keys = df.loc[
                        df[labels[cindex]] == val].groupby(
                        [labels[cindex], labels[output_node_index]]).size()
                    if len(input_nodes_with_prot_keys) > len(
                            output_nodes_with_prot_keys):
                        process_key = '.'.join([val, output_node_val.strip()])
                    elif len(input_nodes_with_prot_keys) < len(
                            output_nodes_with_prot_keys):
                        process_key = '.'.join([input_node_val.strip(), val])
                    else:
                        process_key = '.'.join([input_node_val.strip(), val,
                                                output_node_val.strip()])
                    if process_key not in self.process_map.keys():
                        process = Process(id_=process_key)
                        self.process_map[process_key] = process
                    process_key_sequence.append(process_key)
                elif label.startswith(TableParser.NODE_LABELS):
                    process_key_sequence.append(
                        '.'.join([self._clean_label(label), val]))
            process_key_sequences.append(process_key_sequence)
        for process_key_sequence in process_key_sequences:
            for left, right in self._pairwise(process_key_sequence):
                if left.startswith(TableParser.NODE_LABELS) and not \
                        right.startswith(TableParser.NODE_LABELS):
                    try:
                        material = self.node_map[left]
                    except KeyError:
                        continue
                    process = self.process_map[right]
                    if material not in process.inputs:
                        process.inputs.append(material)
                elif not left.startswith(TableParser.NODE_LABELS) and \
                        right.startswith(TableParser.NODE_LABELS):
                    process = self.process_map[left]
                    try:
                        material = self.node_map[right]
                    except KeyError:
                        continue
                    if material not in process.outputs:
                        process.outputs.append(material)
        for process_key_sequence in process_key_sequences:
            process_only_key_sequence = filter(
                lambda x: not x.startswith(TableParser.NODE_LABELS),
                process_key_sequence)
            for left, right in self._pairwise(process_only_key_sequence):
                left_process = self.process_map[left]
                right_process = self.process_map[right]
                plink(left_process, right_process)

    def _get_value(self, base_column, object_column_group, row):
        cell_value = row[base_column]
        if cell_value == '':
            return cell_value, None
        column_index = list(object_column_group).index(base_column)
        try:
            offset_1r_col = object_column_group[column_index + 1]
            offset_2r_col = object_column_group[column_index + 2]
        except IndexError:
            return cell_value, None
        if offset_1r_col.startswith('Term Source REF') and \
           offset_2r_col.startswith('Term Accession Number'):
            value = OntologyAnnotation(term=str(cell_value))
            term_source_value = row[offset_1r_col]
            if term_source_value is not '':
                try:
                    value.term_source = self.ontology_sources[term_source_value]
                except KeyError:
                    log.debug('term source: %s not found', term_source_value)
            term_accession_value = row[offset_2r_col]
            if term_accession_value is not '':
                value.term_accession = str(term_accession_value)
            return value, None
        try:
            offset_3r_col = object_column_group[column_index + 3]
        except IndexError:
            return cell_value, None
        if offset_1r_col.startswith('Unit') and offset_2r_col.startswith(
                'Term Source REF') \
                and offset_3r_col.startswith('Term Accession Number'):
            category_key = row[offset_1r_col]
            try:
                unit_term_value = self.unit_categories[category_key]
            except KeyError:
                unit_term_value = OntologyAnnotation(term=category_key)
                self.unit_categories[category_key] = unit_term_value
                unit_term_source_value = row[offset_2r_col]
                if unit_term_source_value is not '':
                    try:
                        unit_term_value.term_source = self.ontology_sources[
                            unit_term_source_value]
                    except KeyError:
                        log.debug('term source: %s not found', unit_term_source_value)
                term_accession_value = row[offset_3r_col]
                if term_accession_value is not '':
                    unit_term_value.term_accession = term_accession_value
            return cell_value, unit_term_value
        else:
            return cell_value, None

    def _parse(self, filebuffer):
        raise NotImplementedError(
            'Inherit from this class and implement this method')


class IsaTabSeries(pd.Series):
    @property
    def _constructor(self):
        return IsaTabSeries


class IsaTabDataFrame(pd.DataFrame):

    def __init__(self, *args, **kw):
        super(IsaTabDataFrame, self).__init__(*args, **kw)

    @property
    def _constructor(self):
        return IsaTabDataFrame

    _constructor_sliced = IsaTabSeries

    @staticmethod
    def _clean_label(label):
        for clean_label in TableParser.ALL_LABELS:
            if label.strip().lower().startswith(clean_label.lower()):
                return clean_label
            elif TableParser.RX_CHARACTERISTICS.match(label):
                return 'Characteristics[{val}]'.format(
                    val=next(iter(
                        TableParser.RX_CHARACTERISTICS.findall(label))))
            elif TableParser.RX_PARAMETER_VALUE.match(label):
                return 'Parameter Value[{val}]'.format(
                    val=next(iter(
                        TableParser.RX_PARAMETER_VALUE.findall(label))))
            elif TableParser.RX_FACTOR_VALUE.match(label):
                return 'Factor Value[{val}]'.format(
                    val=next(iter(TableParser.RX_FACTOR_VALUE.findall(label))))
            elif TableParser.RX_COMMENT.match(label):
                return 'Comment[{val}]'.format(
                    val=next(iter(TableParser.RX_COMMENT.findall(label))))

    @property
    def isatab_header(self):
        return list(map(lambda x: self._clean_label(x), self.columns))


class StudySampleTableParser(TableParser):

    def __init__(self, isa=None):
        TableParser.__init__(self)
        if not isinstance(isa, Investigation):
            raise IOError('You must provide an Investigation object output '
                          'from the Investigation parser')
        self.isa = isa
        self.sources = []
        self.samples = []
        self.characteristic_categories = dict()
        self.unit_categories = dict()
        self.process_sequence = []

    def _parse_object_characteristics(self, labels, row):
        characteristics = set()
        for base_column in [x for x in labels if x.startswith(
                ('Characteristics', 'Material Type'))]:
            if base_column.startswith('Material Type'):
                category_key = 'Material Type'
            else:
                category_key = base_column[16:-1]
            try:
                category = self.characteristic_categories[
                    category_key]
            except KeyError:
                category = OntologyAnnotation(term=category_key)
                self.characteristic_categories[category_key] = category
            characteristic = Characteristic(category=category)
            value, unit = self._get_value(base_column, labels, row)
            characteristic.value = value
            characteristic.unit = unit
            characteristics.add(characteristic)
        return sorted(list(characteristics), key=lambda x: x.category.term)

    def _parse_object_comments(self, labels, row):
        comments = set()
        for base_column in (x for x in labels if x.startswith('Comment')):
            comment = Comment(name=base_column[8:-1])
            value, unit = self._get_value(base_column, labels, row)
            comment.value = value
            comment.unit = unit
            comments.add(comment)
        return sorted(list(comments), key=lambda x: x.name)

    def _parse_object_factor_values(self, labels, row):
        factor_values = set()
        for base_column in (x for x in labels if x.startswith('Factor Value')):
            factor_value = FactorValue(factor_name=StudyFactor(
                name=base_column[13:-1]))  # TODO: Check is in Study obj
            value, unit = self._get_value(base_column, labels, row)
            factor_value.value = value
            factor_value.unit = unit
            factor_values.add(factor_value)
        return sorted(list(factor_values), key=lambda x: x.factor_name.name)

    def _parse_materials(self, material_df):
        materials = []
        for _, row in material_df.drop_duplicates().iterrows():
            material_label = next(iter(material_df.columns))
            if material_label.startswith('Source Name'):
                material = Source(name=row[material_label])
            elif material_label.startswith('Sample Name'):
                material = Sample(name=row[material_label])
            material.characteristics = self._parse_object_characteristics(
                material_df.columns, row)
            material.comments = self._parse_object_comments(
                material_df.columns, row)
            material.factor_values = self._parse_object_factor_values(
                material_df.columns, row)
            materials.append(material)
        materials_dict = dict(
            map(lambda x: ('.'.join([material_label, x.name]), x), materials))
        return materials_dict

    def _parse_data_files(self, data_file_df):
        data_files = []
        for _, row in data_file_df.drop_duplicates().iterrows():
            data_file_label = next(iter(data_file_df.columns))
            if data_file_label.startswith('Raw Data File'):
                data_file = RawDataFile(filename=row[data_file_label])
            elif data_file_label.startswith('Derived Data File'):
                data_file = DerivedDataFile(filename=row[data_file_label])
            elif data_file_label.startswith('Derived Spectral Data File'):
                data_file = DerivedSpectralDataFile(
                    filename=row[data_file_label])
            elif data_file_label.startswith('Derived Array Data File'):
                data_file = DerivedArrayDataFile(filename=row[data_file_label])
            elif data_file_label.startswith('Array Data File'):
                data_file = ArrayDataFile(filename=row[data_file_label])
            elif data_file_label.startswith('Protein Assignment File'):
                data_file = ProteinAssignmentFile(filename=row[data_file_label])
            elif data_file_label.startswith('Peptide Assignment File'):
                data_file = PeptideAssignmentFile(filename=row[data_file_label])
            elif data_file_label.startswith(
                    'Post Translational Modification Assignment File'):
                data_file = PostTranslationalModificationAssignmentFile(
                    filename=row[data_file_label])
            elif data_file_label.startswith('Acquisition Parameter Data File'):
                data_file = AcquisitionParameterDataFile(
                    filename=row[data_file_label])
            elif data_file_label.startswith('Free Induction Decay Data File'):
                data_file = FreeInductionDecayDataFile(
                    filename=row[data_file_label])
            # elif data_file_label.startswith('Image File'):
            #     data_file = ImageFile(filename=row[data_file_label])
            # elif data_file_label.startswith('Metabolite Assignment File'):
            #     data_file = MetaboliteAssignmentFile(filename=row[data_file_label])
            elif data_file_label.startswith('Raw Spectral Data File'):
                data_file = RawSpectralDataFile(filename=row[data_file_label])
            data_file.comments = self._parse_object_comments(
                data_file_df.columns, row)
            data_files.append(data_file)
        data_files_dict = dict(
            map(lambda x: ('.'.join([data_file_label, x.name]), x), data_files))
        return data_files_dict

    def _parse_protocol_ref(self, protocol_ref_df):
        pass

    def _parse(self, filebuffer):
        isa_df = IsaTabDataFrame(
            pd.read_csv(filebuffer, dtype=str, sep='\t', encoding='utf-8',
                        comment='#').replace(np.nan, ''))
        isatab_header = isa_df.isatab_header
        object_index = [
            i for i, x in enumerate(isatab_header) if
            x in TableParser.MATERIAL_LABELS + TableParser.DATA_FILE_LABELS +
            tuple(['Protocol REF'])]
        object_columns_grouped = []
        prev_i = next(iter(object_index))
        for curr_i in object_index:
            if prev_i != curr_i:
                object_columns_grouped.append(isa_df.columns[prev_i:curr_i])
            prev_i = curr_i
            object_columns_grouped.append(isa_df.columns[prev_i:])
        for object_column_group in object_columns_grouped:
            object_label = next(iter(object_column_group))  # indicates object
            chopped_isa_df = isa_df[object_column_group]  # chopped is cut by column. sliced is cut by row
            if object_label.startswith('Source Name'):
                sources = self._parse_materials(chopped_isa_df)
                self.node_map.update(sources)
            elif object_label.startswith('Sample Name'):
                samples = self._parse_materials(chopped_isa_df)
                self.node_map.update(samples)
        self.sources = [v for k, v in self.node_map.items() if
                        k.startswith('Source Name')]
        self.samples = [v for k, v in self.node_map.items() if
                        k.startswith('Sample Name')]
        self._make_process_sequence(isa_df)
        self.process_sequence = list(self.process_map.values())


class AssayTableParser(TableParser):

    def __init__(self, isa=None):
        TableParser.__init__(self)
        if not isinstance(isa, Investigation):
            raise IOError('You must provide an Investigation object output '
                          'from the Investigation parser')
        self.isa = isa
        self.samples = None
        self.other_material = None
        self.data_files = None
        self.node_map = dict()
        self.process_map = dict()
        self.process_sequence = None

    def _parse(self, filebuffer):
        df = pd.read_csv(filebuffer, dtype=str, sep='\t', encoding='utf-8',
                         comment='#').replace(np.nan, '')
        samples = dict(
            map(lambda x: ('.'.join(['Sample Name', x]), Sample(name=x)),
                [str(x) for x in df['Sample Name'].drop_duplicates()
                 if x != '']))

        data_files = dict()
        for data_col in (x for x in df.columns if
                         x in TableParser.DATA_FILE_LABELS):
            filenames = [x for x in df[data_col].drop_duplicates() if x != '']
            data_files.update(
                dict(map(
                    lambda x: ('.'.join([data_col, x]),
                               DataFile(filename=x, label=data_col)),
                    filenames)))

        other_material = dict()
        for material_col in (x for x in df.columns if
                             x in TableParser.OTHER_MATERIAL_LABELS):
            if material_col == 'Extract Name':
                extracts = dict(
                    map(lambda x: (
                        '.'.join(['Extract Name', x]), Extract(name=x)),
                        [str(x) for x in df['Extract Name'].drop_duplicates()
                         if x != '']))
                other_material.update(extracts)
            elif material_col == 'Labeled Extract Name':
                labeled_extracts = dict(
                    map(lambda x: (
                        '.'.join(['Labeled Extract Name', x]),
                        LabeledExtract(name=x)),
                        [str(x) for x in
                         df['Labeled Extract Name'].drop_duplicates()
                         if x != '']))
                other_material.update(labeled_extracts)

        self.samples = list(samples.values())
        self.node_map.update(samples)
        self.data_files = list(data_files.values())
        self.node_map.update(data_files)
        self.other_material = list(other_material.values())
        self.node_map.update(other_material)
        self._make_process_sequence(df)
        self.process_sequence = list(self.process_map.values())


class Parser(AbstractParser):

    def __init__(self):
        self.investigation_parser = InvestigationParser()

    def _parse(self, filebuffer):
        self.investigation_parser.parse(filebuffer)

        for study in self.investigation_parser.isa.studies:
            study_sample_table_parser = StudySampleTableParser(
                self.investigation_parser.isa)
            study_sample_table_parser.parse(
                os.path.join(os.path.dirname(filebuffer.name), study.filename))
            study.sources = study_sample_table_parser.sources
            study.samples = study_sample_table_parser.samples
            study.process_sequence = study_sample_table_parser.process_sequence
            for assay in study.assays:
                assay_table_parser = AssayTableParser(
                    self.investigation_parser.isa)
                assay_table_parser.parse(
                    os.path.join(os.path.dirname(filebuffer.name),
                                 assay.filename))
                assay.samples = assay_table_parser.samples
                assay.data_files = assay_table_parser.data_files
                assay.other_material = assay_table_parser.other_material
                assay.process_sequence = assay_table_parser.process_sequence
        self.isa = self.investigation_parser.isa
