"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are characterized by a repeating isoprene unit (C5H8) with a terminal alcohol, phosphate, or diphosphate group.

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

    # Check for terminal group with oxygen or phosphorus
    terminal_group_pattern = Chem.MolFromSmarts("[CX4][OX2,OP(=O)([O-])[O-]]")
    if not mol.HasSubstructMatch(terminal_group_pattern):
        return False, "No terminal oxygen or phosphate/diphosphate group found"
    
     # Check for the presence of the isoprene unit.
    isoprene_unit = Chem.MolFromSmarts("C[C]=[C][C]")
    if not mol.HasSubstructMatch(isoprene_unit):
        return False, "No isoprene unit found"


    # Check if molecule has at least 2 rotatable bonds to verify some kind of chain is present
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
         return False, "Chain too short, less than 2 rotatable bonds"

    
    return True, "Matches prenol criteria"