"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are highly oxygenated tetracyclic triterpenoids with specific structural features.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated core structure for cucurbitacins with potential oxy groups
    # Adjusted to include some typical oxy groups and flexible stereochemistry
    core_pattern = Chem.MolFromSmarts('[C,C@H]1[C,C@H]2[C,C]3[C,C]([C,C]4[C,C@H](=CC(=O)[C,C]4)[C](=O)[C,C](C)[C,C](O)C)C([C,C@H]([C,C]1C)=O)C=C2O') 
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match enhanced tetracyclic skeleton pattern"

    # Check for higher oxygenation
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if oxygen_count / carbon_count < 0.15:  # Ensure sufficient oxygen presence
        return False, "Insufficient oxygenation characteristic of cucurbitacins"

    # Ensure presence of specific double bonds or degree of unsaturation
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count < 3:  # At least 3 double bonds expected 
        return False, "Insufficient unsaturation (double bonds) characteristic of cucurbitacins"

    return True, "Contains enhanced features characteristic of cucurbitacins"