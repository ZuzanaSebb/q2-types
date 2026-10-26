# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import re

from qiime2.plugin import Properties, ValidationError

from q2_types._util import iter_fasta_ids
from q2_types.feature_data import FeatureData
from q2_types.per_sample_sequences import ContigSequencesDirFmt

from .. import Contig
from ...plugin_setup import plugin

_UUID_PATTERN = r"[0-9a-fA-F]{8}-(?:[0-9a-fA-F]{4}-){3}[0-9a-fA-F]{12}"
_SHORTUUID_PATTERN = r"[2-9A-HJ-NP-Za-km-z]{22}"
_PREFIXED_CONTIG_ID = re.compile(
    rf"(?P<sample>.+).(?P<uuid>{_UUID_PATTERN}|{_SHORTUUID_PATTERN})"
)


@plugin.register_validator(FeatureData[Contig % Properties("pooled")])
def _check_pooled_contig_ids(data: ContigSequencesDirFmt, level):
    samples = data.sample_dict()
    if len(samples) != 1:
        raise ValidationError(
            "Pooled contigs must be stored in a single FASTA file; found "
            f"{len(samples)}."
        )

    (filepath,) = samples.values()
    max_records = {"min": 100, "max": None}[level]

    for i, contig_id in enumerate(iter_fasta_ids(filepath)):
        if max_records is not None and i >= max_records:
            break

        if _PREFIXED_CONTIG_ID.fullmatch(contig_id) is None:
            raise ValidationError(
                f'Contig ID "{contig_id}" is not of the form '
                '"<sample-id><separator><uuid or shortuuid>".'
            )
