"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:57925 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion has a protonated amino group ([NH3+]) and a deprotonated carboxyl group ([O-])
    attached to the same alpha-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the zwitterion pattern: [NH3+] attached to a carbon, which is also attached to [O-]
    # This pattern is more flexible to account for cyclic structures and different protonation states
    zwitterion_pattern = Chem.MolFromSmarts("[NH3+,NH2+,NH+][CX4H0,CX3H1][CX3](=[OX1])[OX1-]")
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "No zwitterion pattern found"

    # Verify that the central carbon (alpha-carbon) is attached to both [NH3+] and [O-]
    # This pattern accounts for both linear and cyclic amino acids
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4H0,CX3H1]([NH3+,NH2+,NH+])([CX3](=[OX1])[OX1-])")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "No alpha-carbon with both amino and carboxyl groups found"

    # Check that the alpha-carbon has exactly 4 bonds (including H)
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    for match in alpha_carbon_matches:
        alpha_carbon_idx = match[0]  # Index of the alpha-carbon
        alpha_carbon_atom = mol.GetAtomWithIdx(alpha_carbon_idx)
        
        # Alpha-carbon should have exactly 4 bonds (including H)
        if alpha_carbon_atom.GetTotalDegree() != 4:
            return False, "Alpha-carbon has incorrect number of bonds"

    return True, "Contains alpha-carbon with both amino and carboxyl groups in zwitterionic form"