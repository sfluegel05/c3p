"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:18374 phytosterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    A phytosterol is a plant sterol similar to cholesterol, varying only in carbon side chains
    and/or presence or absence of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general steroid nucleus SMARTS pattern (cyclopentanoperhydrophenanthrene skeleton)
    steroid_nucleus_smarts = '[#6]1([#6])[#6][#6]2[#6]1[#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6][#6]4'
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)

    # Check if molecule contains steroid nucleus
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "Steroid nucleus not found"

    # Identify the steroid skeleton match
    matches = mol.GetSubstructMatches(steroid_nucleus)
    if not matches:
        return False, "Steroid nucleus not found"

    # Check for typical phytosterol side chains at C17 (position after the steroid nucleus)
    # Phytosterols have an alkyl side chain at C17 longer than cholesterol (which has 8 carbons in side chain)
    # We will check for side chains with at least 9 carbons

    # Assuming that the atom at index 17 in the match is C17
    for match in matches:
        steroid_atoms = set(match)
        c17_atom_idx = match[-1]  # Last atom in the steroid nucleus SMARTS
        c17_atom = mol.GetAtomWithIdx(c17_atom_idx)
        # Get side chain atoms connected to C17 that are not part of steroid nucleus
        side_chain_atoms = set()
        atoms_to_visit = []
        for neighbor in c17_atom.GetNeighbors():
            if neighbor.GetIdx() not in steroid_atoms:
                atoms_to_visit.append(neighbor.GetIdx())
                side_chain_atoms.add(neighbor.GetIdx())
        # Traverse the side chain
        while atoms_to_visit:
            current_atom_idx = atoms_to_visit.pop()
            current_atom = mol.GetAtomWithIdx(current_atom_idx)
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in steroid_atoms and neighbor_idx not in side_chain_atoms:
                    side_chain_atoms.add(neighbor_idx)
                    atoms_to_visit.append(neighbor_idx)
        # Count carbon atoms in the side chain
        carbon_count = sum(1 for idx in side_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        # Phytosterols typically have side chains with more than 8 carbons
        if carbon_count >= 9:
            # Optionally check for double bonds in side chain or ring system
            # For simplicity, we will accept this as a phytosterol
            return True, "Molecule contains steroid nucleus with phytosterol characteristic side chain"

    return False, "Side chain at C17 does not match phytosterol characteristics"