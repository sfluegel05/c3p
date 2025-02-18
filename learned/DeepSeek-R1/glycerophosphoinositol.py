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

    # Check for phosphate group (P with double bond to O and three O neighbors)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX2])[OX2][OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check each phosphate for connections to glycerol and inositol
    for p_match in phosphate_matches:
        p_idx = p_match[0]
        p_atom = mol.GetAtomWithIdx(p_idx)

        # Get oxygen atoms connected to phosphorus (excluding double bond)
        o_neighbors = [n for n in p_atom.GetNeighbors() 
                       if n.GetSymbol() == "O" 
                       and mol.GetBondBetweenAtoms(p_idx, n.GetIdx()).GetBondType() != Chem.BondType.DOUBLE]

        glycerol_connected = False
        inositol_connected = False

        for o_atom in o_neighbors:
            # Follow oxygen to connected atom (non-phosphorus)
            connected_atoms = [a for a in o_atom.GetNeighbors() if a.GetIdx() != p_idx]
            if not connected_atoms:
                continue
            target_atom = connected_atoms[0]

            # Check if connected to glycerol backbone (C-O-P)
            # Glycerol pattern: three carbons with two ester groups and one phosphate
            glycerol_candidate = target_atom
            glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-O-C(=O)-*)O-")
            if mol.GetSubstructMatch(glycerol_pattern):
                glycerol_connected = True

            # Check if connected to inositol (cyclohexane with >=4 hydroxyls)
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if len(ring) == 6 and target_atom.GetIdx() in ring:
                    oh_count = 0
                    for atom_idx in ring:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == "O" and atom.GetTotalNumHs() >= 1:
                            oh_count += 1
                    if oh_count >= 4:
                        inositol_connected = True

        if glycerol_connected and inositol_connected:
            # Check for two ester groups (fatty acids)
            ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2][CX3](=O)"))
            if len(ester_matches) >= 2:
                return True, "Contains glycerol backbone with two fatty acid esters, phosphate linked to inositol"

    return False, "Missing required structural components"