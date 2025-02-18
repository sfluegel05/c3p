"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol has a general formula H-[CH2C(Me)=CHCH2]nOH with n > 1 isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one hydroxyl group
    hydroxyls = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1]
    if not hydroxyls:
        return False, "No hydroxyl group found"

    # Check for at least two isoprene units (each with a methyl group and double bond)
    # SMARTS pattern for a carbon with a methyl group and adjacent double bond
    isoprene_carbon = Chem.MolFromSmarts("[C&$(C([CH3])-&!@[C]=,*)]=C")
    matches = mol.GetSubstructMatches(isoprene_carbon)
    if len(matches) < 2:
        return False, f"Found {len(matches)} isoprene units, need at least 2"

    # Verify hydroxyl is attached to the isoprene chain
    for o in hydroxyls:
        carbon = o.GetNeighbors()[0]
        # Check if this carbon is part of a chain with isoprene units
        # Traverse the chain to find connected isoprene carbons
        visited = set()
        stack = [(carbon, 0)]
        max_depth = 0
        while stack:
            atom, depth = stack.pop()
            if atom in visited:
                continue
            visited.add(atom)
            if depth > max_depth:
                max_depth = depth
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor not in visited:
                    stack.append((neighbor, depth + 1))
        # If the chain is long enough (at least 10 atoms for n=2)
        if max_depth >= 9:
            return True, "Contains hydroxyl group and at least two isoprene units in chain"

    return False, "Hydroxyl not part of a sufficient isoprene chain"