"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: CHEBI:7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone has a hydroxy group at the 7-position of the isoflavone core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Corrected SMARTS pattern for 7-hydroxyisoflavone core:
    # Hydroxyl (O) attached to aromatic benzene (A ring) at position 7
    # Fused to chromen-4-one system with substituent (B ring) at position 3
    pattern = Chem.MolFromSmarts("[OH]c1ccc2c(c1)oc(-*)cc(=O)c2")
    if mol.HasSubstructMatch(pattern):
        return True, "7-hydroxy group present on isoflavone core"
    
    # Additional check for cases where the B ring might be connected through different bonds
    # This accounts for possible keto-enol tautomerism or alternate substituent attachments
    alt_pattern = Chem.MolFromSmarts("[OH]c1ccc2c(c1)oc(=O)cc2-*")
    if mol.HasSubstructMatch(alt_pattern):
        return True, "7-hydroxy group present on isoflavone core (alternate pattern)"
    
    return False, "No 7-hydroxy group detected on isoflavone core"