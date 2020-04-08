# Copyright 2020 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from abc import ABC
from collections import defaultdict


class QueryNEW(ABC):
    def __init__(self, label):
        self.__label = label
        self.__gold_standard_df = None
        self.__precision_df = None
        self.__recall_df = None
        self.__confusion_df = None

    @property
    def label(self):
        return self.__label

    @property
    def gold_standard_df(self):
        return self.__gold_standard_df

    @property
    def precision_df(self):
        return self.__precision_df

    @property
    def recall_df(self):
        return self.__recall_df

    @property
    def confusion_df(self):
        return self.__confusion_df

    @label.setter
    def label(self, label):
        self.__label = label

    @gold_standard_df.setter
    def gold_standard_df(self, gold_standard_df):
        self.__gold_standard_df = gold_standard_df

    @precision_df.setter
    def precision_df(self, precision_df):
        self.__precision_df = precision_df

    @recall_df.setter
    def recall_df(self, recall_df):
        self.__recall_df = recall_df

    @confusion_df.setter
    def confusion_df(self, confusion_df):
        self.__confusion_df = confusion_df


class GenomeQueryNEW(QueryNEW):
    binning_type = 'genome'

    def __init__(self, df, label):
        super().__init__(label)
        self.__df = df

    @property
    def df(self):
        return self.__df

    @df.setter
    def df(self, df):
        self.__df = df


class TaxonomicQueryNEW(QueryNEW):
    tax_id_to_parent = None
    tax_id_to_rank = None
    tax_id_to_name = None
    tax_id_to_tax_id = None
    binning_type = 'taxonomic'

    def __init__(self, rank_to_df, label):
        super().__init__(label)
        self.__rank_to_df = rank_to_df

    @property
    def rank_to_df(self):
        return self.__rank_to_df

    @rank_to_df.setter
    def rank_to_df(self, rank_to_df):
        self.__rank_to_df = rank_to_df
