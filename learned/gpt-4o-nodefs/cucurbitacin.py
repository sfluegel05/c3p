"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are highly oxygenated tetracyclic triterpenes with specific structural features.
    
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

    # Define broader patterns for cucurbitacin core (tetracyclic triterpene backbone)
    # Here we don't specify stereochemistry or saturation at connection points
    tetracyclic_pattern = Chem.MolFromSmarts('C1CCC2(CC1)C3CCC4(C2=CC3)C=CC(=O)C4')  # Broad pattern without stereo info
    if not mol.HasSubstructMatch(tetracyclic_pattern):
        return False, "Does not match broad tetracyclic skeleton pattern"

    # Check for high oxygenation
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if oxygen_count / carbon_count < 0.15:  # Adjusted threshold based on known cucurbitacin ratios
        return False, "Insufficient oxygenation characteristic of cucurbitacins"

    # Ensure presence of multiple double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count < 3:  # Ensuring the presence of unsaturation typical of cucurbitacins
        return False, "Insufficient unsaturation (double bonds) characteristic of cucurbitacins"

    # Further features and checks may be added for higher specificity
    return True, "Contains broad features characteristic of cucurbitacins"