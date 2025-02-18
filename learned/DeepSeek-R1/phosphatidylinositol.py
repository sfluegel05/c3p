"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: CHEBI:18179 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol consists of inositol connected via a phosphate group to a
    glycerol backbone which has two fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for inositol (six-membered carbon ring with at least five hydroxyls)
    inositol_found = False
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) != 6:
            continue
        # Check if all atoms in the ring are carbons
        all_carbon = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)
        if not all_carbon:
            continue
        # Count hydroxyl groups attached to the ring
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    oh_count += 1
                    break  # Count each OH once per ring atom
        if oh_count >= 5:
            inositol_found = True
            break
    if not inositol_found:
        return False, "No inositol ring with at least five hydroxyls found"

    # Check for phosphate group connected to inositol and glycerol
    phosphate_pattern = Chem.MolFromSmarts('[O]P(=O)([O])[O]')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group present"

    # Check for glycerol backbone with two ester groups and phosphate connected via oxygen
    # Updated SMARTS pattern to include oxygen between glycerol and phosphate
    glycerol_pattern = Chem.MolFromSmarts('[CH2]OC(=O)-[CH](-OC(=O)-)-[CH2]O-P(=O)(O)')
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with two esters and phosphate not found"

    # Verify that the phosphate connects to both inositol and glycerol
    # Find all atoms in the glycerol pattern
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    for match in glycerol_matches:
        # The phosphate atom is the last atom in the match (index 5)
        p_idx = match[5]
        p_atom = mol.GetAtomWithIdx(p_idx)
        # Check if this phosphate is connected to inositol
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for bond in neighbor.GetBonds():
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetAtomicNum() == 6 and other_atom.IsInRing():
                        # Check if this is part of the inositol ring
                        if any(other_atom.GetIdx() in ring for ring in rings if len(ring) ==6):
                            return True, "Phosphatidyl group esterified to inositol via phosphate linkage"

    return False, "Phosphate not connecting inositol to glycerol"