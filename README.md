# sourmash_plugin_hashannot: extract and/or annotate sourmash matches

[sourmash](https://sourmash.bio) is a tool for biological sequence analysis and comparisons.

The hashannot plugin provides functionality to do two commonly requested
sourmash tasks:

1. extract genomic regions around matching sourmash matches
2. annotate hashes based on genome annotations (GFF files)

## Installation

```
pip install sourmash_plugin_xyz
```

## Usage

The hashannot plugin provides the following sourmash commands (under
`sourmash scripts`):

1. `extract_surrounding` - use sourmash hashes to extract genomic regions from a FASTA file.
2. ...

### `extract_surrounding`

The `extract_surrounding` command can be used to get the underlying
sequence that drives matches between genomes.  It extracts the maximal
contiguous extent of sequence around each matching hash, stopping at
(and not including) neighboring non-matching hashes.

To use it, create a sketch from one genome or set of matches, and then
use that sketch on other genomes.

The following example will extract matching regions from shew-63, using
hashes from shew-47, and place the regions in `shew-matches.fa.gz`.
```
sourmash sketch dna tests/test-data/shew-47.500k.fa.gz -o shew-47.sig.zip
sourmash scripts extract_surrounding tests/test-data/shew-63.500k.fa.gz \
    shew-47.sig.zip -o shew-matches.fa.gz
```

CTB: Future improvements will include greedy merging of regions as
long as the overall match stays above a given containment threshold.

## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention!

## Dev docs

`xyz` is developed at https://github.com/sourmash-bio/sourmash_plugin_template.

### Testing

Run:
```
pytest tests
```

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
make dist
```

followed by `twine upload dist/...`.
