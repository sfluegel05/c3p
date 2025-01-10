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

    # Check for phenolic groups but be flexible with count
    phenolic_pattern = Chem.MolFromSmarts('c:c:c:c:c:c')
    phenolic_matches = mol.GetSubstructMatches(phenolic_pattern)
    if len(phenolic_matches) < 2:
        return False, f"Found only {len(phenolic_matches)} phenolic-like structures, need at least 2"

    # Check for ester and glycosidic-like linkages, broaden patterns
    linkage_patterns = [
        Chem.MolFromSmarts('C(=O)O'),  # ester bond
        Chem.MolFromSmarts('O[C@H]1[C@H](O)C[C@H](O)[C@@H](O)[C@@H]1O'),  # glycosidic linkage
        Chem.MolFromSmarts('c-c(=O)-o'),  # approximation of ester in a ring
    ]

    # Flag for found linkages
    linkage_found = any(mol.HasSubstructMatch(pattern) for pattern in linkage_patterns)
    
    if not linkage_found:
        return False, "Lacking typical tannin linkages"

    return True, "Molecule matches tannin profile with multiple phenolic-like structures and appropriate linkages"