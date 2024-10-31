from rdkit import Chem
from rdkit.Chem import AllChem

def is_sultone(smiles: str):
    """
    Determines if a molecule is a sultone (cyclic ester of hydroxy sulfonic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sultone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for sulfur atom with double bonded oxygens (SO2)
    sulfur_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            # Check if sulfur has double bonded oxygens
            double_bonded_o = 0
            single_bonded_o = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_bonded_o += 1
                    elif bond.GetBondType() == Chem.BondType.SINGLE:
                        single_bonded_o += 1
            if double_bonded_o == 2 and single_bonded_o >= 1:
                sulfur_atoms.append(atom)

    if not sulfur_atoms:
        return False, "No SO2 group found"

    # Check if SO2 is part of a ring
    rings = mol.GetRingInfo()
    
    for sulfur in sulfur_atoms:
        for ring in rings.AtomRings():
            if sulfur.GetIdx() in ring:
                # Check if ring contains O-S bond
                for ring_atom_idx in ring:
                    atom = mol.GetAtomWithIdx(ring_atom_idx)
                    if atom.GetSymbol() == 'O':
                        if mol.GetBondBetweenAtoms(sulfur.GetIdx(), atom.GetIdx()):
                            ring_size = len(ring)
                            return True, f"Found sultone ring of size {ring_size}"

    return False, "SO2 group not part of cyclic ester structure"
# Pr=1.0
# Recall=1.0