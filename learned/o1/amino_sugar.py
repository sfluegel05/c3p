"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:17754 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "Molecule does not contain any rings"

    # Collect candidate sugar rings (5 or 6 membered rings with one oxygen)
    candidate_rings = []
    for ring in rings:
        # Check ring size (5 or 6 members)
        if len(ring) == 5 or len(ring) == 6:
            # Count oxygen atoms in the ring
            num_oxygen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if num_oxygen == 1:
                # Ensure other atoms are carbons
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() in [6, 8] for idx in ring):
                    candidate_rings.append(ring)

    if not candidate_rings:
        return False, "Does not contain a sugar ring"

    # For each candidate ring, check for amino groups replacing hydroxyl groups
    for ring in candidate_rings:
        amino_group_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                # Check neighbors for nitrogen atoms
                neighbors = atom.GetNeighbors()
                for neighbor in neighbors:
                    if neighbor.GetAtomicNum() == 7:  # Nitrogen atom
                        # Check if bond is single (to exclude amides, etc.)
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            amino_group_found = True
                            break
                if amino_group_found:
                    break
        if amino_group_found:
            return True, "Contains sugar ring with amino group(s) replacing hydroxy group(s)"

    return False, "Does not contain amino sugar ring"