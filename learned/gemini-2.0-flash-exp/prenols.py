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

    # Check for at least one terminal alcohol, phosphate or diphosphate group
    terminal_group_pattern = Chem.MolFromSmarts("[OX2H,OP(=O)([O-])[O-]]")
    if not mol.HasSubstructMatch(terminal_group_pattern):
        return False, "No terminal alcohol, phosphate or diphosphate group found"

    # Check for the presence of the basic isoprene unit.
    isoprene_unit = Chem.MolFromSmarts("C[CX3](C)=[CX3][CX4]")
    if not mol.HasSubstructMatch(isoprene_unit):
      return False, "No isoprene unit found"
    
    # Count carbons, hydrogens, and oxygens:
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for sufficient carbons for the isoprene unit(s). Prenols have chains of isoprenes, so the number of carbon atoms should be a multiple of 5 (isoprene) plus some carbons in the terminal group
    if carbon_count < 5 :
      return False, "Too few carbons for an isoprene unit"


    # Check chain length based on rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Chain too short, less than 2 rotatable bonds"
    
    # Approximate number of isoprene units
    isoprene_units = (carbon_count) / 5  
    if isoprene_units < 0.8:
        return False, f"Too few isoprene units ({isoprene_units})."

    return True, "Matches prenol criteria"