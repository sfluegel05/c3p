"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem
from rdkit.Chem import rdchem

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
        if len(ring) == 6:
            ketone_found = False
            double_bond_count = 0

            # Iterate over all bonds in the ring
            for bond_idx in ring:
                bond = mol.GetBondWithIdx(bond_idx)
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    begin_atom_num = bond.GetBeginAtom().GetAtomicNum()
                    end_atom_num = bond.GetEndAtom().GetAtomicNum()
                    
                    # Count any bond as part of a double bond (C=C or C=O)
                    if {begin_atom_num, end_atom_num} == {6, 8}:
                        ketone_found = True
                    else:
                        double_bond_count += 1
            
            # Correctly classify based on presence of ketone and at least one other double bond
            if ketone_found and double_bond_count >= 1:
                return True, "Contains a six-membered ring with a ketone group and at least one other double bond"
            
    return False, "Does not meet cyclohexenone criteria"