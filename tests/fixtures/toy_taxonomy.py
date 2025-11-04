#!/usr/bin/env python3
"""
Toy taxonomy fixture for testing.

Creates a minimal but realistic taxonomy tree for integration tests.
This allows testing reclassify_read() and other functions that require
a real taxonomy tree without needing the full NCBI taxonomy.

Toy Taxonomy Structure:
=======================
1 (root) - no rank
└── 2 (Bacteria) - superkingdom
    └── 1224 (Proteobacteria) - phylum
        └── 561 (Escherichia) - genus
            └── 562 (Escherichia coli) - species

This represents a realistic lineage path from species to root,
allowing us to test:
- Tree traversal
- Confidence calculation at different levels
- LCA resolution
- Lineage caching
"""

from pathlib import Path
from stringmeup.taxonomy import TaxonomyTree


def create_toy_taxonomy_files(output_dir: Path) -> tuple[Path, Path]:
    """
    Create minimal names.dmp and nodes.dmp files for testing.

    Args:
        output_dir: Directory to write the taxonomy files

    Returns:
        Tuple of (names_file_path, nodes_file_path)
    """
    names_file = output_dir / "names.dmp"
    nodes_file = output_dir / "nodes.dmp"

    # Write names.dmp
    # Format: tax_id | name | unique_name | name_class |
    names_content = """1\t|\troot\t|\t\t|\tscientific name\t|
2\t|\tBacteria\t|\tBacteria <bacteria>\t|\tscientific name\t|
1224\t|\tProteobacteria\t|\t\t|\tscientific name\t|
561\t|\tEscherichia\t|\t\t|\tscientific name\t|
562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|
"""
    names_file.write_text(names_content)

    # Write nodes.dmp
    # Format: tax_id | parent_tax_id | rank | ...other_fields... |
    # Using minimal format that TaxonomyTree expects
    nodes_content = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
2\t|\t1\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
1224\t|\t2\t|\tphylum\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
561\t|\t1224\t|\tgenus\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
562\t|\t561\t|\tspecies\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
"""
    nodes_file.write_text(nodes_content)

    return names_file, nodes_file


def create_toy_taxonomy(output_dir: Path) -> TaxonomyTree:
    """
    Create and load a toy taxonomy tree for testing.

    Args:
        output_dir: Directory to write temporary taxonomy files

    Returns:
        TaxonomyTree instance loaded with toy taxonomy

    Example:
        >>> import tempfile
        >>> from pathlib import Path
        >>> tmp = Path(tempfile.mkdtemp())
        >>> taxonomy = create_toy_taxonomy(tmp)
        >>> taxonomy.taxonomy[562].name
        'Escherichia coli'
    """
    names_file, nodes_file = create_toy_taxonomy_files(output_dir)
    return TaxonomyTree(str(nodes_file), str(names_file))


# Taxonomy reference for tests
TOY_TAXONOMY_INFO = {
    'nodes': {
        1: {'name': 'root', 'rank': 'no rank', 'parent': 1},
        2: {'name': 'Bacteria', 'rank': 'superkingdom', 'parent': 1},
        1224: {'name': 'Proteobacteria', 'rank': 'phylum', 'parent': 2},
        561: {'name': 'Escherichia', 'rank': 'genus', 'parent': 1224},
        562: {'name': 'Escherichia coli', 'rank': 'species', 'parent': 561},
    },
    'lineages': {
        562: [562, 561, 1224, 2, 1],  # E. coli to root
        561: [561, 1224, 2, 1],        # Escherichia to root
        1224: [1224, 2, 1],            # Proteobacteria to root
        2: [2, 1],                      # Bacteria to root
        1: [1],                         # Root
    }
}


def get_expected_lineage(tax_id: int) -> list[int]:
    """
    Get the expected lineage for a tax_id in the toy taxonomy.

    Args:
        tax_id: Taxonomic ID

    Returns:
        List of tax IDs from tax_id to root

    Example:
        >>> get_expected_lineage(562)
        [562, 561, 1224, 2, 1]
    """
    return TOY_TAXONOMY_INFO['lineages'].get(tax_id, [])


def get_expected_parent(tax_id: int) -> int:
    """
    Get the expected parent tax_id in the toy taxonomy.

    Args:
        tax_id: Taxonomic ID

    Returns:
        Parent tax ID

    Example:
        >>> get_expected_parent(562)
        561
    """
    return TOY_TAXONOMY_INFO['nodes'].get(tax_id, {}).get('parent')
