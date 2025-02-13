"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:26712 indole alkaloid
An alkaloid containing an indole skeleton
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is an alkaloid containing an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for indole substructure
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)cnc2")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole substructure found"

    # Look for basic nitrogen atom (alkaloid criteria)
    basic_nitrogen_pattern = Chem.MolFromSmarts("[N;+;!X3]")
    if not mol.HasSubstructMatch(basic_nitrogen_pattern):
        return False, "No basic nitrogen atom found (alkaloid criteria)"

    # Check molecular weight - alkaloids typically <700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 700:
        return False, "Molecular weight too high for alkaloid"

    # Count nitrogen atoms - alkaloids typically have 1-3
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1 or n_count > 3:
        return False, "Number of nitrogen atoms outside expected range for alkaloids"

    # Count aromatic rings - alkaloids typically have 2-4
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 2 or aromatic_rings > 4:
        return False, "Number of aromatic rings outside expected range for alkaloids"

    return True, "Contains an indole substructure and meets alkaloid criteria"