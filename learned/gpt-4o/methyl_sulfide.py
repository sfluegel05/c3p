"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide where at least one of the groups attached to sulfur is a methyl group (CH3-S).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the methyl sulfide pattern (sulfur bonded directly to a methyl group)
    methyl_group = Chem.MolFromSmarts("[CH3]")
    sulfur_atom = Chem.MolFromSmarts("[S]")
    
    for atom in mol.GetAtoms():
        # Check if the atom is a sulfur atom
        if atom.GetSymbol() == 'S':
            for neighbor in atom.GetNeighbors():
                # Check if any neighbor is a methyl group
                if neighbor.GetSymbol() == 'C' and mol.HasSubstructMatch(methyl_group):
                    if neighbor.GetTotalNumHs() == 3:
                        return True, "Contains a sulfur atom bonded to a methyl group"
    
    return False, "Does not contain a sulfur atom bonded to a methyl group"

# Example usage:
# smiles = "CCSC"  # Example SMILES for ethyl methyl sulfide
# result, reason = is_methyl_sulfide(smiles)
# print(result, reason)