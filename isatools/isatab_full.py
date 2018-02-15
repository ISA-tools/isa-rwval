
from __future__ import absolute_import

import os

import numpy as np
import pandas as pd

from isatools.isatab_meta import AbstractParser, InvestigationParser, TableParser
from isatools.model import *


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
