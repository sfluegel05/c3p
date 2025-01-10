"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH
    in which the carbon skeleton is composed of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; prenols are acyclic"
    
    # Check for hydroxyl group (-OH) attached to a carbon chain
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No terminal hydroxyl group (-OH) found attached to carbon"
    
    # Ensure only allowed atoms are present (C, H, O)
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Atom {atom.GetSymbol()} not allowed in prenols"
    
    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons % 5 != 0:
        return False, f"Number of carbon atoms ({num_carbons}) is not a multiple of 5"
    
    # Check for repeating isoprene units
    # Define a pattern for the isoprene unit: C=C-C-C
    # This pattern allows for variations in double bond positions due to stereochemistry
    isoprene_unit_pattern = Chem.MolFromSmarts("C=C-C-C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit_pattern)
    num_isoprene_units = len(isoprene_matches)
    if num_isoprene_units == 0:
        return False, "No isoprene units found"
    
    # Verify that the number of isoprene units corresponds to the number of carbons
    expected_isoprene_units = num_carbons // 5
    if num_isoprene_units < expected_isoprene_units:
        return False, f"Expected at least {expected_isoprene_units} isoprene units, found {num_isoprene_units}"
    
    # All checks passed; classify as prenol
    return True, f"Molecule is a prenol with {num_isoprene_units} isoprene units and a terminal hydroxyl group"

__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'prenols',
        'definition': "Any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of one or more isoprene units (biogenetic precursors of the isoprenoids).",
        'parents': []
    }
}