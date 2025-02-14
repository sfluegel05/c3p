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

    # Check for isoprene unit [C]([C])=[C][C] with methyl attached
    isoprene_unit = Chem.MolFromSmarts("[C]([C])=[C][C]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit)
    if len(isoprene_matches) < 1:
        return False, "No isoprene unit found"

    # Check for alcohol terminal group
    alcohol_terminal = Chem.MolFromSmarts("[OX2H1]")
    alcohol_match = mol.HasSubstructMatch(alcohol_terminal)

    # Check for phosphate terminal group
    phosphate_terminal = Chem.MolFromSmarts("[P](=[O])([O-])([O-])")
    phosphate_match = mol.HasSubstructMatch(phosphate_terminal)

    # Check for diphosphate terminal group
    diphosphate_terminal = Chem.MolFromSmarts("[P]([O])(=[O])([O][P]([O])(=[O])[O-])")
    diphosphate_match = mol.HasSubstructMatch(diphosphate_terminal)

    # Check if at least one terminal group is present
    if not (alcohol_match or phosphate_match or diphosphate_match):
       return False, "No terminal alcohol, phosphate or diphosphate group found"
    
    # Check for the presence of at least two isoprene units (for prenol to be "repeating")
    if len(isoprene_matches) < 2:
         return False, "Less than 2 isoprene units"
    
    # Check if molecule has at least 2 rotatable bonds to verify some kind of chain is present
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
         return False, "Chain too short, less than 2 rotatable bonds"

    return True, "Matches prenol criteria"