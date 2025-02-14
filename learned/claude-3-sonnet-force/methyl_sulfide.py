"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:33023 methyl sulfide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.

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

    # Look for sulfur atom
    sulfur_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atom found"

    # Check if any sulfur atom has a neighboring methyl group
    for sulfur_idx in sulfur_atoms:
        for neighbor_idx in mol.GetAtomWithIdx(sulfur_idx).GetNeighbors():
            neighbor = mol.GetAtomWithIdx(neighbor_idx)
            if neighbor.GetSymbol() == 'C' and len(neighbor.GetNeighbors()) == 4:
                hydrogen_count = sum(1 for neighbor_neighbor_idx in neighbor.GetNeighbors() if mol.GetAtomWithIdx(neighbor_neighbor_idx).GetSymbol() == 'H')
                if hydrogen_count == 3:
                    return True, "Contains a sulfur atom with a neighboring methyl group"

    return False, "No methyl group found attached to sulfur"