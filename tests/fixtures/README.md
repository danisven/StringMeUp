# Test Fixtures for StringMeUp

This directory contains test fixtures for integration and end-to-end tests.

## Toy Taxonomy

### `toy_taxonomy.py`

Helper module for creating a minimal taxonomy tree for testing.

**Taxonomy Structure:**
```
1 (root) - no rank
└── 2 (Bacteria) - superkingdom
    └── 1224 (Proteobacteria) - phylum
        └── 561 (Escherichia) - genus
            └── 562 (Escherichia coli) - species
```

**Usage:**
```python
from tests.fixtures.toy_taxonomy import create_toy_taxonomy
from pathlib import Path

# In a test
def test_with_taxonomy(tmp_path):
    taxonomy = create_toy_taxonomy(tmp_path)
    # Use taxonomy in tests
    assert taxonomy.taxonomy[562].name == 'Escherichia coli'
```

**Functions:**
- `create_toy_taxonomy(output_dir)`: Create and load a TaxonomyTree
- `create_toy_taxonomy_files(output_dir)`: Create names.dmp and nodes.dmp files
- `get_expected_lineage(tax_id)`: Get expected lineage for a tax ID
- `get_expected_parent(tax_id)`: Get expected parent for a tax ID

## Example Kraken2 Files

### `example.kraken2`

Standard single-end Kraken2 output with 10 reads covering various scenarios:

- **read1**: 60% hits, 30% unclassified, 10% ambiguous (threshold crossing scenario)
- **read2**: 80% hits, 20% unclassified (high confidence)
- **read3**: 40% hits, 60% unclassified (low confidence)
- **read4**: Unclassified read (100% unclassified)
- **read5**: Hits to genus level (561) + some species (562)
- **read6**: Hits across multiple levels (species, genus, phylum)
- **read7**: 90% hits, 10% ambiguous
- **read8**: Hits to superkingdom level (Bacteria)
- **read9**: Perfect classification (100% hits)
- **read10**: 25% hits, 50% unclassified, 25% ambiguous

### `example_paired.kraken2`

Paired-end Kraken2 output with 5 read pairs:

- **read1**: Paired reads with ambiguous in both
- **read2**: Paired reads, no ambiguous
- **read3**: Mixed genus/species hits across pair
- **read4**: Unclassified pair
- **read5**: High ambiguous content in both reads

**Format:** Uses `|:|` delimiter between paired reads

### `example_with_ambiguous.kraken2`

Focused on ambiguous k-mer scenarios (v0.2.0 bug fix verification):

- **read_10pct_ambiguous**: 10% ambiguous k-mers
- **read_20pct_ambiguous**: 20% ambiguous k-mers
- **read_50pct_ambiguous**: 50% ambiguous k-mers (extreme case)
- **read_only_ambiguous**: 100% ambiguous (edge case)
- **read_no_ambiguous**: 0% ambiguous (baseline)
- **read_threshold_cross**: Example of threshold crossing due to ambiguous

## Expected Outputs

### Confidence Calculations with Toy Taxonomy

Using the toy taxonomy and confidence threshold of 0.5:

| Read | Original | K-mer String | Confidence at 562 | Expected Reclassification |
|------|----------|--------------|-------------------|---------------------------|
| read1 | 562 (E. coli) | 562:60 0:30 A:10 | 0.60 (60/100) | 562 (stays) |
| read2 | 562 (E. coli) | 562:80 0:20 | 0.80 (80/100) | 562 (stays) |
| read3 | 562 (E. coli) | 562:40 0:60 | 0.40 (40/100) | Moves up (genus/phylum/bacteria) |
| read5 | 561 (Escherichia) | 561:50 562:20 0:30 | 0.70 (70/100) | 561 (stays) |
| read6 | 562 (E. coli) | 562:30 561:20 1224:10 0:40 | At 562: 0.30, At 561: 0.50, At 1224: 0.60 | 1224 or higher |

### Critical v0.2.0 Bug Fix Example

**Scenario**: read_threshold_cross with threshold 0.65

- K-mer string: `562:60 0:30 A:10`
- **Old buggy calculation**: 60 / 90 = 0.667 → PASSES threshold
- **New correct calculation**: 60 / 100 = 0.60 → FAILS threshold
- **Result**: Read must move up the tree (correct behavior)

## Usage in Tests

### Example Integration Test

```python
import pytest
from pathlib import Path
from tests.fixtures.toy_taxonomy import create_toy_taxonomy
from stringmeup.stringmeup import reclassify_read, ReadClassification

def test_reclassification_with_toy_data(tmp_path):
    """Test reclassification using toy taxonomy and example data."""
    # Load toy taxonomy
    taxonomy = create_toy_taxonomy(tmp_path)

    # Create read from example data
    read = ReadClassification(
        id="read1",
        original_taxid=562,
        kmer_string="562:60 0:30 A:10"
    )

    # Reclassify
    result, _ = reclassify_read(
        read=read,
        confidence_threshold=0.5,
        taxonomy_tree=taxonomy,
        verbose_input=False,
        minimum_hit_groups=1,
        taxa_lineages={},
        paired_input=False
    )

    # Verify
    assert result.classified == True
    assert result.ambiguous_kmers == 10  # v0.2.0 bug fix
    # Confidence = 60 / 100 = 0.60 (includes ambiguous)
```

## Adding New Fixtures

When adding new test fixtures:

1. **Name descriptively**: `example_<scenario>.kraken2`
2. **Document in this README**: Add description and expected behavior
3. **Keep realistic**: Base on actual Kraken2 output format
4. **Cover edge cases**: Include unusual but valid scenarios
5. **Provide expected outputs**: Document what should happen

## Maintenance

- Fixtures should be **minimal** but **representative**
- Update when Kraken2 format changes
- Keep aligned with toy taxonomy structure
- Version control all fixture files
