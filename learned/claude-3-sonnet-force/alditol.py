"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:17937 alditol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is a carbohydrate that is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for acyclic structure
    if mol.HasSubstructMatch(Chem.MolFromSmarts("R1R2R3R4")):
        return False, "Molecule contains rings, but alditols are acyclic"
    
    # Check for polyol structure (multiple hydroxy groups)
    num_hydroxy = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and atom.GetHybridization() == Chem.HybridizationType.SP3)
    if num_hydroxy < 3:
        return False, "Fewer than 3 hydroxy groups, not a polyol"
    
    # Check for terminal hydroxy groups
    terminal_hydroxy = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetDegree() == 1)
    if terminal_hydroxy != 2:
        return False, "Does not have 2 terminal hydroxy groups as required for alditols"
    
    # Count carbon atoms
    num_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if number of hydroxy groups matches the expected formula (n + 2)
    if num_hydroxy != num_carbon + 2:
        return False, f"Incorrect number of hydroxy groups for alditol (expected {num_carbon + 2}, got {num_hydroxy})"

    return True, "Acyclic polyol matching the formula HOCH2[CH(OH)]nCH2OH"