"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    Polyamines typically contain multiple amine groups which may be dispersed 
    in various structural contexts across the molecule.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure presence of multiple amine groups
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;+0,!$([N]C=O)]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)

    if len(amine_matches) < 2:
        return False, f"Found {len(amine_matches)} amine group(s), need at least 2"

    # Check that the structure supports a polyamine nature (flexible interpretation)
    if len(amine_matches) >= 3 or mol.GetRingInfo().NumRings() > 0:
        return True, "Contains multiple amine groups supporting a polyamine structure"

    # Consider small non-linear or interconnected structural motifs
    significant_connections = [
        (i, j) for i in range(len(amine_matches))
        for j in range(i + 1, len(amine_matches))
        if Chem.rdmolops.GetShortestPath(mol, amine_matches[i][0], amine_matches[j][0])  # significant path
    ]

    if significant_connections:
        return True, "Contains spaced amine groups or belongs to common polyamine motifs"

    return False, "Amine groups are not sufficiently matched to polyamine characteristics"