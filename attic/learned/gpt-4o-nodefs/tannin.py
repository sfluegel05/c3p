"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are large polyphenolic compounds composed of phenolic acids
    or catechins, often forming complex structures through ester or
    glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for multiple aromatic rings (common in tannins: catechin or ellagic units)
    aromatic_ring_pattern = Chem.MolFromSmarts('c1ccccc1')
    aromatic_matches = mol.GetSubstructMatches(aromatic_ring_pattern)
    if len(aromatic_matches) < 2:
        return False, f"Found only {len(aromatic_matches)} aromatic rings, need at least 2"

    # Check for a variety of linking bond types to assert tannin-like complexity
    linkage_patterns = [
        Chem.MolFromSmarts('C(=O)O'),                              # Ester bond
        Chem.MolFromSmarts('c1cc(O)c(O)cc1'),                      # phenolic rings with hydroxyl groups
        Chem.MolFromSmarts('O[C@H]1[C@H](O)C[C@H](O)[C@@H](O)[C@@H]1O'),  # Glycosidic linkage
        Chem.MolFromSmarts('c-c(=O)-o'),                          # Approximation of ester in a ring
    ]

    # Flag for recognized linkages
    linkage_found = any(mol.HasSubstructMatch(pattern) for pattern in linkage_patterns)
    
    if not linkage_found:
        return False, "Lacking typical tannin linkages or catechin hinges"

    return True, "Molecule matches tannin profile with multiple aromatic rings and appropriate linkages"