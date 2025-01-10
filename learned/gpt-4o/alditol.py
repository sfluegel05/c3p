"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2-[CH(OH)]n-CH2OH, 
    derived from an aldose by reduction of the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule is acyclic
    if rdmolops.GetSSSR(mol) > 0:
        return False, "Molecule contains a ring structure, not an acyclic alditol"

    # Define the SMARTS pattern for the alditol structure
    alditol_pattern = Chem.MolFromSmarts("OCC(O)[CH2]")  # Main unit in repeating sequence
    if alditol_pattern is None:
        return (None, "Invalid SMARTS pattern")

    # Ensure molecule ends with 'CH2OH'
    end_group_pattern = Chem.MolFromSmarts("CO")
    start_group_pattern = Chem.MolFromSmarts("CO[CH2]")

    # Check full structure for alditol pattern matching
    match = mol.GetSubstructMatches(alditol_pattern)
    if len(match) < 2:  # Must be at least two repeating C(OH) units between CH2OH ends
        return False, f"Structure with {len(match)} C-OH units, insufficient for alditol"

    # Verify the ends of the molecule have the appropriate groups
    if not mol.HasSubstructMatch(end_group_pattern):
        return False, "Molecule does not end with a 'CH2OH' group"
    if not mol.HasSubstructMatch(start_group_pattern):
        return False, "Molecule does not begin with a 'HOCH2' group"

    return True, "SMILES string matches the structural pattern of an acyclic alditol"