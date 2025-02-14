"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    for ring in ring_info.BondRings():
        if len(ring) == 6:  # Focus on six-membered rings
            ketone_count = 0
            other_double_bond_count = 0

            # Iterate over all bonds in the ring
            for bond_idx in ring:
                bond = mol.GetBondWithIdx(bond_idx)
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    begin_atom = bond.GetBeginAtom()
                    end_atom = bond.GetEndAtom()

                    # Check if the bond is part of a ketone (C=O)
                    if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or \
                       (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6):
                        ketone_count += 1
                    else:
                        # Count other types of double bonds
                        other_double_bond_count += 1
            
            # Verify that there is exactly one ketone and exactly one additional double bond
            if ketone_count == 1 and other_double_bond_count == 1:
                return True, "Contains a six-membered ring with one ketone group and another double bond"
            
    return False, "Does not meet cyclohexenone criteria"