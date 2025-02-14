"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: CHEBI:27370 gamma-lactone

A gamma-lactone is a lactone with a five-membered lactone ring.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find lactone rings
    lactone_pattern = Chem.MolFromSmarts("[O;r5]1[C;r5][C;r5][C;r5][C;r5]1=O")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    
    # Check if any match is a 5-membered lactone ring
    has_gamma_lactone = False
    for match in lactone_matches:
        ring = Chem.rdmolops.GetMolFragFromAtomSmilesPattern(mol, smiles, atomsToUse=match, bondsToUse=[], atomMemory={})
        if ring is not None and ring.GetNumAtoms() == 5:
            has_gamma_lactone = True
            break
    
    if has_gamma_lactone:
        return True, "Contains a five-membered lactone ring"
    else:
        return False, "No five-membered lactone ring found"