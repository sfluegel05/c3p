"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:34084 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a (potential) aldehyde group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 5 or 6 carbon atoms (allowing for some flexibility)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5 or c_count > 6:
        return False, f"Expected 5 or 6 carbons, found {c_count}"

    # Check for aldehyde group (either explicit or potential)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    # If no explicit aldehyde, check for potential aldehyde (e.g., in cyclic form)
    if not aldehyde_matches:
        # Look for a carbon with a single bond to oxygen (potential aldehyde in cyclic form)
        potential_aldehyde_pattern = Chem.MolFromSmarts("[CX4][OX2]")
        potential_aldehyde_matches = mol.GetSubstructMatches(potential_aldehyde_pattern)
        if not potential_aldehyde_matches:
            return False, "No aldehyde or potential aldehyde group found"

    # Check for multiple hydroxyl groups (at least 3) attached to the carbon backbone
    hydroxyl_pattern = Chem.MolFromSmarts("[C][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Expected at least 3 hydroxyl groups attached to the carbon backbone, found {len(hydroxyl_matches)}"

    # Check molecular weight (should be around 150 g/mol for aldopentoses)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 200:
        return False, f"Molecular weight {mol_wt:.2f} is outside expected range for aldopentoses"

    return True, "Contains 5 or 6 carbons, a (potential) aldehyde group, and multiple hydroxyl groups attached to the carbon backbone"