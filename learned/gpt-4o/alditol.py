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

    # Define the SMARTS pattern for an alditol structure
    # The pattern should match [CH2OH][CH(OH)]n[CH2OH] with n >= 2
    alditol_pattern = Chem.MolFromSmarts("OCC(O)C(O)*CO")
    if alditol_pattern is None:
        return (None, "Invalid SMARTS pattern")

    # Check if the molecule contains cycles
    if rdmolops.GetSSSR(mol) > 0:
        return False, "Molecule contains a ring structure, not an acyclic alditol"

    # Check if the molecule matches the alditol SMARTS pattern
    match = mol.HasSubstructMatch(alditol_pattern)
    if not match:
        return False, "No continuous acyclic polyol chain matching alditol found"

    return True, "SMILES string matches the structural pattern of an alditol"