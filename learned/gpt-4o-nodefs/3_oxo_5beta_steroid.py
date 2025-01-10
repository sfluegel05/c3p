"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid should have a steroid backbone with a ketone group at the 3-position
    and specific stereochemistry denoted as '5beta'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for 3-oxo group (simple ketone pattern, needs contextualization)
    oxo_pattern = Chem.MolFromSmarts("C(=O)C")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not found"

    # Basic pattern for a steroid backbone (simplified; needs actual steroid recognition)
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C)CCC4C(C)CCCC4=C3C2=C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid backbone not found"

    # Check for chiral centers and configuration at position 5 (beta orientation)
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    has_5beta = any(center[0] == 5 and '@' in center[1] for center in chiral_centers)
    if not has_5beta:
        return False, "5beta stereochemistry not confirmed; chiral centers: {}".format(chiral_centers)

    return True, "Molecule matches the 3-oxo-5beta-steroid characteristics"