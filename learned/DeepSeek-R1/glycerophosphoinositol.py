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

    # Check for phosphate group (at least two oxygen substituents)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check glycerol backbone: three carbons with two ester groups and one phosphate
    # Glycerol pattern: C1-C2-C3 where C3 is connected to phosphate via O
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-[OX2]C(=O))-[CH2]-[OX2]-[PX4]")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
        return False, "Glycerol backbone with two esters and phosphate not found"

    # Verify phosphate is connected to inositol (cyclohexane with >=4 hydroxyls)
    inositol_found = False
    for p_match in phosphate_matches:
        p_atom = mol.GetAtomWithIdx(p_match[0])
        # Check oxygen neighbors of phosphate for connection to inositol
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O':
                # Follow O to check if part of inositol ring
                next_atoms = [a for a in neighbor.GetNeighbors() if a.GetIdx() != p_match[0]]
                if not next_atoms:
                    continue
                next_atom = next_atoms[0]
                # Check if part of a 6-membered ring with multiple OHs
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

    # Check for exactly two ester groups (from fatty acids) on glycerol
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2]C(=O)"))
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    return True, "Contains glycerol with two fatty acid esters, phosphate linked to inositol"