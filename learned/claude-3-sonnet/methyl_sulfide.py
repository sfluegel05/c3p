"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:33108 methyl sulfide
Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide with at least one methyl group attached to the sulfur.

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
    
    # Find sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    
    # Check if any sulfur atom has a methyl group attached
    has_methyl_sulfide = False
    for sulfur in sulfur_atoms:
        for neighbor in sulfur.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if carbon is part of a methyl group (has 3 hydrogen neighbors and is not part of a ring)
                if len([n for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 1]) == 3 and not neighbor.IsInRing():
                    has_methyl_sulfide = True
                    break
        if has_methyl_sulfide:
            break
    
    if not has_methyl_sulfide:
        return False, "No methyl group attached to sulfur"
    
    # Check if the molecule is an aliphatic sulfide
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic atoms, not an aliphatic sulfide"
    
    if len(sulfur_atoms) == 0:
        return False, "Does not contain a sulfur atom, not a sulfide"
    
    return True, "Contains aliphatic sulfide with methyl group attached to sulfur"