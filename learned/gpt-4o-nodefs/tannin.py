"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds typically featuring multiple phenolic groups.

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

    # Look for aromatic rings with at least 2 hydroxyl groups
    aromatic_phenol_pattern = Chem.MolFromSmarts('c1cc([OH])c([OH])ccc1')
    if not mol.HasSubstructMatch(aromatic_phenol_pattern):
        return False, "No aromatic phenolic rings with sufficient hydroxyl groups found"
    
    # Look for multiple phenolic rings
    phenolic_ring_count = len(mol.GetSubstructMatches(aromatic_phenol_pattern))
    if phenolic_ring_count < 2:
        return False, f"Found only {phenolic_ring_count} phenolic rings, need at least 2"

    # Look for glycosidic links (sugar-like structures) or ester bonds
    sugar_pattern = Chem.MolFromSmarts('O1CC(O)C(O)C(O)C(O)C1')  # Example SMARTS for a pyran sugar ring
    ester_pattern = Chem.MolFromSmarts('COC(=O)C')
    
    has_glycosidic_linkage = mol.HasSubstructMatch(sugar_pattern)
    has_ester_bond = mol.HasSubstructMatch(ester_pattern)
    
    if not (has_glycosidic_linkage or has_ester_bond):
        return False, "No glycosidic or ester linkages found typical of tannins"

    return True, "Contains multiple phenolic groups and necessary linkages typical of tannins"