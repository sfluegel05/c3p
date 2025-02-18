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
    These molecules have a glycerol backbone with two fatty acids, a phosphate group,
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

    # Check for glycerol backbone with two ester groups and one phosphate
    # Glycerol pattern: C-O-C(=O) (two esters) and C-O-P (phosphate)
    glycerol_pattern = Chem.MolFromSmarts(
        "[CH2]OC(=O)[!O][CH](OC(=O)[!O])[CH2]OP(=O)([OX2])[OX2]"
    )
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with two esters and phosphate not found"

    # Find myo-inositol ring (six-membered carbon ring with multiple hydroxyls)
    inositol_rings = []
    for ring in Chem.GetSymmSSSR(mol):
        if len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in ring):
            inositol_rings.append(ring)
    if not inositol_rings:
        return False, "No inositol ring found"

    # Check for at least one phosphate group attached to inositol ring
    phosphate_pattern = Chem.MolFromSmarts("[C][OX2]P(=O)([OX2])[OX2]")
    phosphate_found = any(mol.HasSubstructMatch(phosphate_pattern) for _ in inositol_rings)
    if not phosphate_found:
        return False, "Inositol lacks phosphate groups"

    # Verify minimum hydroxyl count on inositol (at least 4)
    for ring in inositol_rings:
        hydroxyl_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    hydroxyl_count += 1
        if hydroxyl_count >= 4:
            return True, "Contains glycerol backbone with two esters, phosphate-linked myo-inositol with phosphates"

    return False, "Inositol lacks sufficient hydroxyl groups"