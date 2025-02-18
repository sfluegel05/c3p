"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: CHEBI:82525 hydroxynaphthoquinone
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone (naphthalene with two ketone groups) 
    substituted by at least one hydroxyl group on the naphthalene moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find naphthalene systems (two fused benzene rings)
    naphtha_pattern = Chem.MolFromSmarts("c1ccc2ccccc2c1")  # Naphthalene pattern
    naphtha_matches = mol.GetSubstructMatches(naphtha_pattern)
    if not naphtha_matches:
        return False, "No naphthalene moiety found"

    # Collect all atoms in naphthalene systems
    naphtha_atoms = set()
    for match in naphtha_matches:
        naphtha_atoms.update(match)

    # Check for exactly two ketone groups (C=O) in the naphthalene system
    ketone_count = 0
    for atom_idx in naphtha_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:  # Only check carbon atoms
            continue
        # Look for double bonds to oxygen (C=O)
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    ketone_count += 1
    if ketone_count != 2:
        return False, f"Found {ketone_count} ketone groups in naphthalene (needs exactly 2)"

    # Check for hydroxyl groups attached to naphthalene carbons
    hydroxyl_found = False
    hydroxyl_pattern = Chem.MolFromSmarts("[c]-[OH]")  # Hydroxyl attached to aromatic carbon
    for match in mol.GetSubstructMatches(hydroxyl_pattern):
        carbon_idx = match[0]
        if carbon_idx in naphtha_atoms:
            hydroxyl_found = True
            break
    if not hydroxyl_found:
        return False, "No hydroxyl group attached to naphthalene system"

    return True, "Naphthoquinone with hydroxyl substituent on aromatic ring"