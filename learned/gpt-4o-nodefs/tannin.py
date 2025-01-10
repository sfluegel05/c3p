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

    # More flexible search for aromatic phenolic rings
    phenolic_ring_pattern = Chem.MolFromSmarts('c1c([OH])c(O)c(O)cc1')
    if not mol.HasSubstructMatch(phenolic_ring_pattern):
        return False, "No aromatic phenolic rings found with multiple hydroxyl groups"
    
    # Count phenolic rings
    phenolic_ring_count = len(mol.GetSubstructMatches(phenolic_ring_pattern))
    if phenolic_ring_count < 3:
        return False, f"Found only {phenolic_ring_count} phenolic rings, need at least 3"

    # Explore more diverse connectivity: ester and glycosidic linkages
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    glycosidic_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)C[C@H](O)[C@@H](O)[C@@H]1O')
    
    has_ester_bond = mol.HasSubstructMatch(ester_pattern)
    has_glycosidic_linkage = mol.HasSubstructMatch(glycosidic_pattern)
    
    if not (has_ester_bond or has_glycosidic_linkage):
        return False, "Lacking ester or glycosidic linkages typical of tannins"

    return True, "Molecule matches tannin profile with multiple phenolic rings and appropriate linkages"