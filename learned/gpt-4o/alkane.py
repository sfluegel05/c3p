"""
Classifies: CHEBI:18310 alkane
"""
from rdkit import Chem

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is an acyclic, saturated hydrocarbon with the general formula CnH2n+2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if only C and H atoms are present
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 1):  # Carbon is 6, Hydrogen is 1
            return False, "Contains atoms other than C and H"

    # Check that all carbon atoms are sp3 hybridized
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() != Chem.HybridizationType.SP3:
            return False, "Contains unsaturated carbon atoms"

    # Check if the molecule is acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structure"

    # Check CnH2n+2 formula
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if h_count != 2 * c_count + 2:
        return False, f"Formula mismatch: expected H = {2 * c_count + 2}, found H = {h_count}"

    return True, "Valid alkane: meets all criteria for an acyclic, saturated hydrocarbon"