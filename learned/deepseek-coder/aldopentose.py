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

    # Check for exactly 5 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected exactly 5 carbons, found {c_count}"

    # Check for explicit aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No explicit aldehyde group found"

    # Check for exactly 4 hydroxyl groups attached to the carbon backbone
    hydroxyl_pattern = Chem.MolFromSmarts("[C][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) != 4:
        return False, f"Expected exactly 4 hydroxyl groups attached to the carbon backbone, found {len(hydroxyl_matches)}"

    # Check molecular weight (should be around 150 g/mol for aldopentoses)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 140 or mol_wt > 160:
        return False, f"Molecular weight {mol_wt:.2f} is outside expected range for aldopentoses"

    return True, "Contains exactly 5 carbons, an explicit aldehyde group, and exactly 4 hydroxyl groups attached to the carbon backbone"