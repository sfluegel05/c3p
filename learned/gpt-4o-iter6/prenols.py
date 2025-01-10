"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols with the formula H-[CH2C(Me)=CHCH2]nOH, containing one or more isoprene units.

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
    
    # Pattern for isoprene unit with optional methyl branching
    isoprene_pattern = Chem.MolFromSmarts("[CH2]-C(=C)-C[CH3]")
    
    # Search for isoprene units
    isoprene_blocks = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_blocks) < 1:
        return False, "No or insufficient isoprene units found"
    
    # Check for precisely one alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) != 1:
        return False, "There must be exactly one alcohol group"
    
    # Validate the position of OH to be at the terminus of the molecule
    terminal_oxygen = [a.GetIdx() for a in mol.GetAtomsWithQuery(Chem.MolFromSmarts("[OX2H]")) if a.GetDegree() == 1]
    if not terminal_oxygen:
        return False, "Alcohol group must be terminal"

    # Ensure alcohol group is part of the main structure, not isolated
    for atom_idx in terminal_oxygen:
        neighboring_atoms = [n.GetAtomicNum() for n in mol.GetAtomWithIdx(atom_idx).GetNeighbors()]
        if 6 in neighboring_atoms:  # Carbon is present
            return True, "Contains isoprene units with a terminal alcohol group"
    
    return False, "Can't confirm isoprene-alcohol linkage"

__metadata__ = {'chemical_class': {'name': 'prenol'}}