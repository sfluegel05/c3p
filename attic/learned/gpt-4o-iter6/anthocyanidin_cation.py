"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are oxygenated derivatives of flavylium (2-phenylchromenylium).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for characteristic flavylium cation pattern (C1c2cc[cH+]c(O)c2cc(=O)c1) with oxygens
    flavylium_pattern = Chem.MolFromSmarts("[cH+]1cc2ccccc2c(c(=O)c1)O")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "Does not match the flavylium cation pattern"
    
    # Check for multiple oxygen groups, typically at least 2 or more hydroxyls or ethers
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() > 1)
    if o_count < 3:
        return False, "Insufficient number of oxygen substituents"

    return True, "Contains 2-phenylchromenylium core with sufficient oxygenation"