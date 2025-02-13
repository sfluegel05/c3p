"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:33108 methyl sulfide
Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide with at least one methyl group attached to the sulfur,
    either directly or through a chain of aliphatic carbon atoms.

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
    
    # Check if any sulfur atom has a methyl group attached, either directly or through a chain of aliphatic carbon atoms
    has_methyl_sulfide = False
    for sulfur in sulfur_atoms:
        if has_methyl_group_attached(mol, sulfur):
            has_methyl_sulfide = True
            break
    
    if not has_methyl_sulfide:
        return False, "No methyl group attached to sulfur (either directly or through aliphatic chain)"
    
    # Check if the molecule is an aliphatic sulfide
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic atoms, not an aliphatic sulfide"
    
    if len(sulfur_atoms) == 0:
        return False, "Does not contain a sulfur atom, not a sulfide"
    
    return True, "Contains aliphatic sulfide with methyl group attached to sulfur (either directly or through aliphatic chain)"

def has_methyl_group_attached(mol, sulfur_atom):
    """
    Checks if a sulfur atom has a methyl group attached, either directly or through a chain of aliphatic carbon atoms.
    """
    visited = set()
    queue = [(sulfur_atom, 0)]
    
    while queue:
        atom, depth = queue.pop(0)
        
        if atom in visited:
            continue
        visited.add(atom)
        
        # Check if atom is a carbon with 3 hydrogen neighbors and not in a ring (i.e., methyl group)
        if atom.GetAtomicNum() == 6 and len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]) == 3 and not atom.IsInRing():
            return True
        
        # Check if atom is an aliphatic carbon (not in a ring and not attached to heteroatoms other than hydrogen)
        if atom.GetAtomicNum() == 6 and not atom.IsInRing() and all(n.GetAtomicNum() in (1, 6) for n in atom.GetNeighbors()):
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    queue.append((neighbor, depth + 1))
        
        # Stop searching if depth exceeds a reasonable limit (e.g., 10)
        if depth > 10:
            break
    
    return False