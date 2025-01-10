"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    Proanthocyanidins are characterized by flavan-3-ol units linked together, often with hydroxyphenyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavan-3-ol pattern: C1(C(C(O)C2)C3=c4ccc(O)cc4Oc4ccc(O)cc34)C1
    # Check for aromatic rings with hydroxyl groups
    aromatic_oh_pattern = Chem.MolFromSmarts('c1cc(O)ccc1')
    if not mol.HasSubstructMatch(aromatic_oh_pattern):
        return False, "No aromatic rings with hydroxyl groups found"
    
    # Check if there are more than one linked flavonoid units (e.g., catechin/epicatechin)
    flavonoid_pattern = Chem.MolFromSmarts('C1(C(C(O)C2)C3)C2')
    if len(mol.GetSubstructMatches(flavonoid_pattern)) < 2:
        return False, "Insufficient flavan-3-ol units linked together"
    
    # Ensure multiple chiral centers common in proanthocyanidins
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if num_chiral_centers < 2:
        return False, "Not enough chiral centers"

    return True, "Contains features consistent with a proanthocyanidin: aromatic hydroxyl rings and multiple flavonoid units linked"