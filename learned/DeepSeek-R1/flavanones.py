"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:flavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones have a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define the core dihydrobenzopyran-4-one structure (chromone with dihydro)
    core_pattern = Chem.MolFromSmarts('c12c(C(=O)CC(O1)C)cccc2')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No dihydrobenzopyran-4-one core found"

    # Check for aryl group attached to the C2 position of the core
    # The aryl group is a benzene ring (possibly substituted) attached to the carbon adjacent to the oxygen
    for match in mol.GetSubstructMatches(core_pattern):
        # Find the oxygen atom in the core
        oxygen_atom = None
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8 and atom.IsInRing():
                oxygen_atom = atom
                break
        if not oxygen_atom:
            continue

        # Get the carbon adjacent to the oxygen in the core (C2)
        c_adjacent = None
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetIdx() in match:  # Ensure it's part of the core
                c_adjacent = neighbor
                break
        if not c_adjacent:
            continue

        # Check substituents of C2 for aromatic rings
        for neighbor in c_adjacent.GetNeighbors():
            if neighbor.GetIdx() in match:
                continue  # Skip core atoms
            # Check if the substituent is part of an aromatic ring
            if neighbor.GetIsAromatic():
                return True, "Contains dihydrobenzopyran-4-one core with aryl substituent"
            # Check if the substituent is connected to an aromatic ring
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    return True, "Contains dihydrobenzopyran-4-one core with aryl substituent"
            # Check rings for aromaticity
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if neighbor.GetIdx() in ring:
                    is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                    if is_aromatic and len(ring) >= 6:
                        return True, "Contains dihydrobenzopyran-4-one core with aryl substituent"

    return False, "No aryl substituent found on dihydrobenzopyran-4-one core"