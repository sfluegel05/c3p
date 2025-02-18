"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: CHEBI:xxxxx polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate has a phosphate group attached to the terminal allylic carbon of a polyprenol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define phosphate ester pattern: C-O-P(=O)(O)O
    phosphate_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-P(=O)([OX2-])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester group detected"

    # Check if phosphate is attached to an allylic carbon (adjacent to double bond)
    matches = mol.GetSubstructMatches(phosphate_pattern)
    for match in matches:
        # Carbon attached to phosphate oxygen is the first atom in the match
        carbon_idx = match[0]
        carbon = mol.GetAtomWithIdx(carbon_idx)
        
        # Check neighbors of this carbon for double bonds
        for neighbor in carbon.GetNeighbors():
            # Check all bonds of the neighbor
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Ensure the double bond isn't with the original carbon (conjugated case)
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetIdx() != carbon_idx:
                        return True, "Phosphate attached to allylic carbon in polyprenol chain"

    return False, "Phosphate group not attached to allylic carbon"