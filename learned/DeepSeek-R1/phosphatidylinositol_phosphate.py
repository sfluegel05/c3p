"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: phosphatidylinositol phosphate (PIP)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES.
    These molecules have a glycerol backbone with two fatty acid esters, a phosphate group,
    and a myo-inositol ring with at least one phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PIP, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Improved glycerol pattern: three carbons with two esters and one phosphate
    # Allows any chain configuration and stereochemistry
    glycerol_pattern = Chem.MolFromSmarts(
        "[CH2]OC(=O)-[CH](-OC(=O)-)-[CH2]OP(=O)(O)"
    )
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with two esters and phosphate not found"

    # Find phosphate connected to glycerol
    phosphate_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CH2]OP(=O)(O)"))
    if not phosphate_matches:
        return False, "No phosphate group on glycerol"

    # Track inositol connection through phosphate oxygen
    for p_match in phosphate_matches:
        p_atom = mol.GetAtomWithIdx(p_match[2])  # P atom in CH2-OP group
        # Find oxygen connecting to inositol
        for o_atom in p_atom.GetNeighbors():
            if o_atom.GetAtomicNum() == 8 and o_atom.GetDegree() == 2:
                # Follow connection to inositol ring
                inositol_atom = next((a for a in o_atom.GetNeighbors() if a.GetIdx() != p_match[2]), None)
                if inositol_atom and inositol_atom.GetAtomicNum() == 6:
                    # Check if part of a six-membered carbon ring
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if len(ring) == 6 and inositol_atom.GetIdx() in ring:
                            # Verify inositol has at least one phosphate group (not the connecting one)
                            for atom_idx in ring:
                                atom = mol.GetAtomWithIdx(atom_idx)
                                for bond in atom.GetBonds():
                                    if bond.GetBondType() == Chem.BondType.SINGLE:
                                        neighbor = bond.GetOtherAtom(atom)
                                        if neighbor.GetAtomicNum() == 15:  # Phosphorus
                                            # Check for phosphate group (P=O)
                                            if any(b.GetBondType() == Chem.BondType.DOUBLE for b in neighbor.GetBonds()):
                                                return True, "Contains glycerol backbone with two esters, phosphate-linked myo-inositol with phosphates"

    return False, "No phosphorylated myo-inositol ring found"