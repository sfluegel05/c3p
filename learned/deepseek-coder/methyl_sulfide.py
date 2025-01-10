"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:48375 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is any aliphatic sulfide where at least one of the organyl groups attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl sulfide pattern: sulfur bonded to a methyl group
    # Ensure the sulfur is not in a ring or aromatic
    methyl_sulfide_pattern = Chem.MolFromSmarts("[SX2;!R;!a][CX4;H3]")

    # Check if the molecule contains the methyl sulfide pattern
    matches = mol.GetSubstructMatches(methyl_sulfide_pattern)
    if not matches:
        return False, "No sulfur atom bonded to a methyl group found"

    # Additional checks to ensure proper classification
    for match in matches:
        sulfur_idx = match[0]  # Index of the sulfur atom in the match
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # Check if sulfur is part of an aliphatic chain
        if not sulfur_atom.IsInRing() and not sulfur_atom.GetIsAromatic():
            # Check that the methyl group is directly bonded to sulfur
            methyl_idx = match[1]
            methyl_atom = mol.GetAtomWithIdx(methyl_idx)
            if methyl_atom.GetSymbol() == 'C' and methyl_atom.GetTotalNumHs() == 3:
                # Ensure that the sulfur is only bonded to aliphatic carbons or hydrogen
                for neighbor in sulfur_atom.GetNeighbors():
                    if neighbor.GetSymbol() != 'C' and neighbor.GetSymbol() != 'H':
                        return False, "Sulfur is bonded to non-aliphatic or non-hydrogen atoms"
                return True, "Contains an aliphatic sulfur atom bonded to at least one methyl group"

    return False, "No valid aliphatic methyl sulfide group found"