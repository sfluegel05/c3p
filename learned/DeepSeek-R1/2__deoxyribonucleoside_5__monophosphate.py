"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    Criteria: Deoxyribose (no 2'-OH), phosphate at 5' position, nucleobase at 1' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphate group in 5' position using SMARTS
    # Pattern: deoxyribose (no 2'-OH) with 5'-phosphate and nucleobase
    # [O] is ring oxygen, adjacent C (C2') has no OH, C5' has CH2OP group
    deoxy_pattern = Chem.MolFromSmarts(
        "[O;R1][C@H]1[C@@H]([C@H](O)[C@@H](COP(=O)([OH,O-])[OH,O-])O1)N"
    )
    if mol.HasSubstructMatch(deoxy_pattern):
        return True, "Has deoxyribose with 5'-phosphate and nucleobase"

    # Alternative pattern accounting for different base attachment points
    alt_pattern = Chem.MolFromSmarts(
        "[O;R1][C@@H]1[C@H](O)[C@@H](COP(=O)([OH,O-])[OH,O-])O[C@H]1[NH0]"
    )
    if mol.HasSubstructMatch(alt_pattern):
        return True, "Matches alternative pattern with nucleobase"

    # Check for absence of 2'-OH specifically
    # Get all potential 2' positions (adjacent to ring oxygen in furanose)
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 5:  # Furanose ring check
            continue
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:  # Find ring oxygen
                neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in ring]
                if len(neighbors) != 2:
                    continue
                # Check adjacent carbons (potential C1' and C2')
                for neighbor_idx in neighbors:
                    neighbor = mol.GetAtomWithIdx(neighbor_idx)
                    if neighbor.GetAtomicNum() == 6:
                        # Check if this is C2' (no OH)
                        has_oh = False
                        for bond in neighbor.GetBonds():
                            other_atom = bond.GetOtherAtomIdx(neighbor_idx)
                            if mol.GetAtomWithIdx(other_atom).GetAtomicNum() == 8 and \
                               mol.GetAtomWithIdx(other_atom).GetTotalNumHs() > 0:
                                has_oh = True
                                break
                        if not has_oh:
                            # Check if connected to phosphate and nucleobase
                            phosphate_present = any(
                                a.GetAtomicNum() == 15 for a in neighbor.GetNeighbors()
                            )
                            base_present = any(
                                a.GetAtomicNum() == 7 and a.IsInRing() for a in neighbor.GetNeighbors()
                            )
                            if phosphate_present and base_present:
                                return True, "Structural features match"
    
    return False, "Does not meet all criteria"