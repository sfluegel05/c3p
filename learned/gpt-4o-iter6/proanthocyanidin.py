"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer consisting of two or more hydroxyflavan units, often linked by C4-C8 or C4-C6 bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for phenolic groups in proanthocyanidins
    phenol_pattern = Chem.MolFromSmarts("c1c(O)ccc(O)c1")  # Generalized phenolic pattern
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) < 4:  # Checking for at least 2 phenolic groups per flavan unit, suggests at least 2 flavan units present
        return False, "Insufficient phenolic groups to suggest multiple hydroxyflavan units"
    
    # More generalized C4-C8 and C4-C6 linkage pattern search
    linkage_pattern = Chem.MolFromSmarts("[C@H]1(c2cc(O)cc(O)c2)C1-[c]3[c]([OH])[c]([OH])[c]([OH])[c]([c]3[OH])")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)
    
    if len(linkage_matches) < 1:  # Needs at least one C4-CX linkage
        return False, "No interflavan linkage found"
    
    return True, "Contains multiple hydroxyflavan units with at least one interflavan linkage"