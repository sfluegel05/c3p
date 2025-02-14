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

    # Define the flavone core SMARTS pattern
    # This pattern represents a chromen-4-one core with a 2-aryl substitution
    flavone_smarts = '''
    [$([cH]1[cH][cH][cH][cH][cH]1)]    # Benzene ring fused to pyran ring
    -c2cc(=O)                          # Connection to pyran ring with ketone at position 4
    oc3ccccc23                         # Pyran ring fused to another benzene ring (chromen-4-one)
    '''
    flavone_smarts = flavone_smarts.replace('\n', '').replace(' ', '')
    flavone_mol = Chem.MolFromSmarts(flavone_smarts)

    if flavone_mol is None:
        return False, "Invalid SMARTS pattern for flavone core"

    # Check for flavone core match
    if not mol.HasSubstructMatch(flavone_mol):
        return False, "Does not contain flavone core structure"

    # Now, ensure that the aryl group at position 2 is present
    # Get the matches of the flavone core
    matches = mol.GetSubstructMatches(flavone_mol)
    for match in matches:
        # In the SMARTS pattern, the atom index of position 2 (the carbon connected to the aryl group) is 6
        position2_atom_idx = match[6]
        position2_atom = mol.GetAtomWithIdx(position2_atom_idx)

        # Get neighbors of position 2 atom not in the flavone core
        flavone_core_atom_idxs = set(match)
        neighbors = [nbr for nbr in position2_atom.GetNeighbors() if nbr.GetIdx() not in flavone_core_atom_idxs]

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