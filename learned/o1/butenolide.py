"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the generalized butenolide SMARTS pattern
    butenolide_smarts = '[O]=C1[C;R][C;R]=[C;R][O;R]1'
    pattern = Chem.MolFromSmarts(butenolide_smarts)

    # Check if molecule contains the butenolide core
    if mol.HasSubstructMatch(pattern):
        return True, "Contains butenolide ring (2-furanone skeleton)"
    else:
        # Additional check using ring analysis
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) != 5:
                continue  # Skip if not a five-membered ring

            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Check for one oxygen atom in the ring
            o_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
            if len(o_atoms) != 1:
                continue  # Need exactly one oxygen atom in the ring

            # Check for one carbonyl group (C=O) within the ring
            c_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 6]
            carbonyl_carbons = []
            for atom in c_atoms:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in ring:
                        continue
                    if neighbor.GetAtomicNum() == 8 and atom.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        carbonyl_carbons.append(atom)
            if len(carbonyl_carbons) != 1:
                continue  # Need exactly one carbonyl group in the ring

            # Check that the oxygen atom is adjacent to the carbonyl carbon (forming the lactone)
            lactone_formed = False
            carbonyl_carbon = carbonyl_carbons[0]
            oxygen_atom = o_atoms[0]
            if mol.GetBondBetweenAtoms(carbonyl_carbon.GetIdx(), oxygen_atom.GetIdx()):
                lactone_formed = True

            if not lactone_formed:
                continue

            # Check for at least one double bond in the ring (excluding the carbonyl)
            double_bond_found = False
            for bond in mol.GetBonds():
                if bond.IsInRing() and bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        # Exclude the carbonyl double bond
                        begin_atom = bond.GetBeginAtom()
                        end_atom = bond.GetEndAtom()
                        if not ((begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6)):
                            double_bond_found = True
                            break
            if not double_bond_found:
                continue  # Need at least one double bond in the ring besides the carbonyl

            # All conditions met, it's a butenolide
            return True, "Contains butenolide ring (2-furanone skeleton)"

        return False, "Does not contain butenolide ring (2-furanone skeleton)"