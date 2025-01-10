"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as containing a lipid portion connected via glycosidic linkage to a sugar moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Look for a lipid backbone pattern, generalized to include sphingo-like lipids
    lipid_backbone_patterns = [
        Chem.MolFromSmarts("C(COC(=O)[CX4])"),  # General ester linkage (e.g., glycerolipid)
        Chem.MolFromSmarts("N[C@@H](C)[CX4]"),  # Amide linkage (e.g., typical in sphingolipids)
    ]
    
    if not any(mol.HasSubstructMatch(pat) for pat in lipid_backbone_patterns):
        return False, "No recognizable lipid backbone structure found"

    # Look for glycosidic linkage pattern
    glycosidic_linkage_pattern = Chem.MolFromSmarts("C-O-C")
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No glycosidic linkage detected"

    # Recognize a broad pattern for sugar moieties
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1"),  # Simple monosaccharide unit
        Chem.MolFromSmarts("C1OC(C(O)C(O1))"),  # 1,2 linked sugar ring
    ]

    if not any(mol.HasSubstructMatch(pat) for pat in sugar_patterns):
        return False, "No recognizable sugar moiety found"
    
    return True, "The structure matches the criteria for a glycolipid, containing a lipid backbone linked via glycosidic linkage to a sugar moiety"