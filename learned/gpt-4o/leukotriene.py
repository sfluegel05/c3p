"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is characterized as a C20 polyunsaturated fatty acid with specific conjugated double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the backbone has roughly 20 carbon atoms
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_chain_length < 20:
        return False, "Molecule does not have adequate icosanoid carbon chain length"
    
    # Check for the presence of four double bonds which is essential for leukotrienes
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_matches < 4:
        return False, f"Found {double_bond_matches} double bonds, need at least 4"
    
    # Check for at least three conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    has_conjugation = mol.HasSubstructMatch(conjugated_pattern)
    if not has_conjugation:
        return False, "No conjugated double bonds pattern found"
    
    # Additional group checks, like hydroxyl groups or other typical leukotriene modifications
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    
    # Deduce classification with collected reasons
    classification_reasons = [
        "Has adequate carbon chain length",
        "Has 4 double bonds, of which pattern supports conjugation" if has_conjugation else "Double bond conjugation pattern insufficient",
        f"Contains {oh_matches} hydroxyl groups" if oh_matches > 0 else "No hydroxyl groups detected"
    ]
    reason = "; ".join(filter(None, classification_reasons))
    
    return True, reason

# Example usage:
# is_leukotriene("CCCCCCCCCCCCCCCCCCCCC/C=C/C=C/C=C/C=C/CCCC(=O)O")