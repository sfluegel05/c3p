"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone is defined as a compound containing a 2-aryl-1-benzopyran-4-one
    (2-arylchromen-4-one) skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the chromone core (1-benzopyran-4-one)
    chromone_smarts = 'c1cc2oc(=O)cc2cc1'  # Chromone core pattern
    chromone_mol = Chem.MolFromSmarts(chromone_smarts)

    if chromone_mol is None:
        return False, "Invalid SMARTS pattern for chromone core"

    # Find matches for the chromone core
    matches = mol.GetSubstructMatches(chromone_mol)
    if not matches:
        return False, "Does not contain chromone core"

    # For each match, check if carbon at position 2 is attached to an aryl group
    for match in matches:
        chromone_atoms = match
        # Atom at position 2 is the third atom in the chromone core pattern
        position2_atom_idx = chromone_atoms[2]
        position2_atom = mol.GetAtomWithIdx(position2_atom_idx)

        # Get neighbors of position 2 atom not in chromone core
        neighbors = [nbr for nbr in position2_atom.GetNeighbors() if nbr.GetIdx() not in chromone_atoms]

        # Check if any neighbor is part of an aromatic ring (aryl group)
        for nbr in neighbors:
            if nbr.GetIsAromatic():
                # Check if neighbor is part of an aromatic ring
                ring_info = mol.GetRingInfo()
                atom_rings = ring_info.AtomRings()
                is_aryl = False
                for ring in atom_rings:
                    if nbr.GetIdx() in ring and len(ring) >= 5:
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                            is_aryl = True
                            break
                if is_aryl:
                    return True, "Contains flavone core structure with 2-aryl group"

    return False, "Does not contain flavone core structure with 2-aryl group"