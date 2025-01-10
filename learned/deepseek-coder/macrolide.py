"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"

    # Find largest ring size and its atoms
    atom_rings = ring_info.AtomRings()
    largest_ring = max(atom_rings, key=len)
    max_ring_size = len(largest_ring)
    if max_ring_size < 12:
        return False, f"Largest ring has {max_ring_size} members, need at least 12"

    # Look for lactone group (cyclic ester) in largest ring
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[OX2][#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check if any ester is in the largest ring
    ester_in_largest_ring = any(
        any(atom_idx in largest_ring for atom_idx in match)
        for match in ester_matches
    )
    if not ester_in_largest_ring:
        return False, "No lactone group in largest ring"

    # Check molecular weight (macrolides are typically >400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight {mol_wt:.1f} too low for macrolide"

    # Check for typical macrolide features
    # Count oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, "Too few oxygen atoms for macrolide"

    # Count carbonyl groups (C=O)
    carbonyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and 
                        any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
                            bond.GetOtherAtom(atom).GetAtomicNum() == 8 
                            for bond in atom.GetBonds()))
    
    # Count hydroxyl groups (OH)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and 
                        any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE and 
                            bond.GetOtherAtom(atom).GetAtomicNum() == 1 
                            for bond in atom.GetBonds()))
    
    # Check for polyketide-like structure (alternating carbonyls and hydroxyls)
    # Look for pattern: [C]-[C](=O)-[C]-[C]-OH
    polyketide_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6][#6][#6][OX2H]")
    has_polyketide = mol.HasSubstructMatch(polyketide_pattern)
    
    # Overall oxygen-containing functional groups check
    if carbonyl_count + hydroxyl_count < 2 and not has_polyketide:
        return False, "Insufficient oxygen-containing functional groups for macrolide"

    # Check for typical macrolide carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for typical macrolide"

    # Check molecular complexity
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Molecular structure too simple for macrolide"

    return True, "Contains macrocyclic lactone ring with 12+ members and oxygen-rich structure"