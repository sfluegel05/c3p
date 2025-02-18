"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: CHEBI:XXXXX glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol has a glycerol backbone with two fatty acid esters,
    a phosphate group at sn-3 position, and inositol attached to the phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group connected to glycerol and inositol
    # Phosphate pattern: P(=O)(O)(O)O connected to two different groups
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX2]-[!P])[OX2]-[!P][OX2]")
    if not phosphate_pattern:
        return False, "Invalid phosphate pattern"
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check glycerol backbone: three carbons with two ester groups and one phosphate
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-[OX2]-C(=O))-[CH2]-[OX2]-[P]")
    if not glycerol_pattern:
        return False, "Invalid glycerol pattern"
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
        return False, "Glycerol backbone with two esters and phosphate not found"

    # Verify phosphate is at sn-3 position (third carbon in glycerol)
    # Assuming the glycerol pattern is C1-C2-C3 where C3 is connected to phosphate
    c3_idx = glycerol_match[-3]  # Last three atoms in match are C3, O, P?
    phosphate_atom = mol.GetAtomWithIdx(glycerol_match[-1])
    if phosphate_atom.GetSymbol() != 'P':
        return False, "Phosphate not correctly attached to glycerol"

    # Check inositol (cyclohexane with >=4 hydroxyls) connected to phosphate
    inositol_found = False
    for p_match in phosphate_matches:
        p_idx = p_match[0]
        # Check oxygens connected to P
        for neighbor in mol.GetAtomWithIdx(p_idx).GetNeighbors():
            if neighbor.GetSymbol() == 'O':
                # Follow O to check if connected to inositol ring
                next_atom = [a for a in neighbor.GetNeighbors() if a.GetIdx() != p_idx]
                if not next_atom:
                    continue
                next_atom = next_atom[0]
                # Check if part of a 6-membered ring with multiple OH groups
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if len(ring) == 6 and next_atom.GetIdx() in ring:
                        oh_count = 0
                        for atom_idx in ring:
                            atom = mol.GetAtomWithIdx(atom_idx)
                            if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() >= 1:
                                oh_count += 1
                        if oh_count >= 4:
                            inositol_found = True
                            break
                if inositol_found:
                    break
        if inositol_found:
            break
    if not inositol_found:
        return False, "Inositol ring not found connected to phosphate"

    # Check for exactly two ester groups (from fatty acids)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2][CX3](=O)"))
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    return True, "Contains glycerol with two fatty acid esters, phosphate linked to inositol"