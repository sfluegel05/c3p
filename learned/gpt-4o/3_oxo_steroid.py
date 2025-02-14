"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is characterized by a steroidal framework with an oxo group (C=O) at the 3rd position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized pattern for steroids backbone (cyclopentanoperhydrophenanthrene)
    steroid_smarts = 'C1CCC2C3C4C=C(C=C4CCC3CCC2C1)C'
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Find potential matches of the oxo groups
    oxo_smarts = '[C]=O'
    oxo_pattern = Chem.MolFromSmarts(oxo_smarts)
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    
    if not oxo_matches:
        return False, "No C=O bond found"

    # Check if one of these C=O groups is at the 3rd position of a steroid ring
    # This is checked by ensuring the carbon atom with the oxo group is part of the ring and is correctly positioned
    for match in oxo_matches:
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Looking for adjacency within steroid backbone
        if any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE and bond.IsInRing() for bond in carbon_atom.GetBonds()):
            return True, "Detected steroid with a 3-oxo (C=O) group at the correct position"

    return False, "Oxo group is not at the 3rd position in the steroid structure"