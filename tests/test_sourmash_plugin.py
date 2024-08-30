"""
Tests for sourmash_plugin_xyz.
"""
import os
import pytest

import sourmash
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


def test_run_sourmash(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('', fail_ok=True)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status != 0                    # no args provided, ok ;)


def test_simple(runtmp):
    shew47 = utils.get_test_data('shew-47.500k.fa.gz')
    shew63 = utils.get_test_data('shew-63.500k.fa.gz')

    # sketch 47
    sig47 = runtmp.output('47.sig.zip')
    runtmp.sourmash('sketch', 'dna', shew47, '-o', sig47)

    # extract from 63
    contigs = runtmp.output('out.fa')
    runtmp.sourmash('scripts', 'extract_surrounding', shew63, sig47,
                    '-o', contigs)

    out = runtmp.last_result.out
    assert "original 500000, extracted 293015" in out

    assert os.path.exists(contigs)


def test_simple_identical(runtmp):
    shew47 = utils.get_test_data('shew-47.500k.fa.gz')

    # sketch 47
    sig47 = runtmp.output('47.sig.zip')
    runtmp.sourmash('sketch', 'dna', shew47, '-o', sig47)

    # extract from 47
    contigs = runtmp.output('out.fa')
    runtmp.sourmash('scripts', 'extract_surrounding', shew47, sig47,
                    '-o', contigs)

    out = runtmp.last_result.out
    assert "original 500000, extracted 500000" in out

    assert os.path.exists(contigs)


def test_make_hashannotdb(runtmp):
    shew47 = utils.get_test_data('shew-47.500k.fa.gz')
    shew47_gff = utils.get_test_data('shew-47.500k.gff')
    outdb = runtmp.output('foo.sqldb')

    runtmp.sourmash('scripts', 'make_hashannotdb',
                    shew47, shew47_gff, '-o', outdb)
