"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26964 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is a terpenoid derived from a tetraterpene (C40 skeleton),
    which may be rearranged or modified by the removal of skeletal atoms (generally methyl groups).
    It typically features a long conjugated polyene chain derived from eight isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 50:
        return False, f"Carbon count ({c_count}) not in range for tetraterpenoid (25-50 carbons)"
    
    # Find the longest conjugated polyene chain
    # We consider a chain of alternating single and double bonds between carbons
    mol = Chem.AddHs(mol)
    longest_chain = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            visited = set()
            stack = [(atom, None, 0)]
            while stack:
                current_atom, prev_bond_order, length = stack.pop()
                visited.add(current_atom.GetIdx())
                for bond in current_atom.GetBonds():
                    neighbor = bond.GetOtherAtom(current_atom)
                    if neighbor.GetAtomicNum() != 6:
                        continue
                    if neighbor.GetIdx() in visited:
                        continue
                    bond_order = bond.GetBondTypeAsDouble()
                    if prev_bond_order is None or bond_order != prev_bond_order:
                        stack.append((neighbor, bond_order, length + 1))
                if length > longest_chain:
                    longest_chain = length
    if longest_chain < 12:
        return False, f"Longest conjugated polyene chain length is {longest_chain}, need at least 12"

    # Optional: Check for multiple methyl groups (CH3)
    methyl_group = Chem.MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_group)
    if len(methyl_matches) < 4:
        return False, f"Found {len(methyl_matches)} methyl groups, need at least 4"
    
    # Optional: Approximate number of isoprene units
    # Each isoprene unit is C5, so estimate number of isoprene units
    isoprene_units = c_count / 5
    if isoprene_units < 5 or isoprene_units > 10:
        return False, f"Number of isoprene units ({isoprene_units:.1f}) not in range for tetraterpenoid (5-10 units)"
    
    # Check for terpenoid features: presence of isoprene-like motifs
    isoprene_smarts = "[C](=C)C=C[C]"
    isoprene_mol = Chem.MolFromSmarts(isoprene_smarts)
    isoprene_matches = mol.GetSubstructMatches(isoprene_mol)
    if len(isoprene_matches) < 3:
        return False, f"Found {len(isoprene_matches)} isoprene-like units, need at least 3"
    
    return True, "Molecule meets criteria for tetraterpenoid (derived from tetraterpene with long conjugated polyene chain)"