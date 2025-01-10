"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose must have an aldehyde group and hydroxyl groups typically at every carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aldehyde group - must have terminal C=O
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group (C=O) found"
    
    # Count hydroxyl groups (OH)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if len(oh_matches) < 2:
        return False, "Fewer hydroxyl groups than expected for an aldose"

    # Check for a carbon backbone
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbon atoms for an aldose"

    # Ensure the molecule has chiral centers typical of aldoses like glucose
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 2:
        return False, "Insufficient chiral centers for an aldose"
    
    return True, "Contains aldehyde group with sufficient hydroxyl groups and carbon backbone characteristic of an aldose"