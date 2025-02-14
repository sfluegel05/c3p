"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a glycerol backbone with one fatty acid chain attached via an ester bond and a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphate group (-P(=O)(O)-OH or -P(=O)(O)2 or -P(=O)([O-])(O)) or -P(=O)([O-])([O-])
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O-])([OX2])")
    if not mol.HasSubstructMatch(phosphate_pattern):
       return False, "No phosphate group found"
    
    # Look for 1 ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    #Check that the phosphate group is connected to the glycerol backbone
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    phosphate_match = mol.GetSubstructMatch(phosphate_pattern)
    found_phosphate_link = False
    for g_atom_idx in glycerol_match:
        g_atom = mol.GetAtomWithIdx(g_atom_idx)
        for p_atom_idx in phosphate_match:
            p_atom = mol.GetAtomWithIdx(p_atom_idx)
            for neighbor in g_atom.GetNeighbors():
                if neighbor.GetIdx() == p_atom_idx:
                    found_phosphate_link = True
                    break
        if found_phosphate_link:
            break
    if not found_phosphate_link:
       return False, "Phosphate group not connected to glycerol backbone"
    
     #Check that the ester group is connected to the glycerol backbone
    found_ester_link = False
    for e_match in ester_matches:
        ester_atom_idx = e_match[0]
        ester_atom = mol.GetAtomWithIdx(ester_atom_idx)
        for g_atom_idx in glycerol_match:
            g_atom = mol.GetAtomWithIdx(g_atom_idx)
            for neighbor in g_atom.GetNeighbors():
                if neighbor.GetIdx() == ester_atom_idx:
                    found_ester_link = True
                    break
            if found_ester_link:
                break

    if not found_ester_link:
        return False, "Ester group not connected to glycerol backbone"

    
     # Check that the phosphate and the ester are attached to different carbons on the glycerol backbone
    glycerol_phosphate_atom = None
    for g_atom_idx in glycerol_match:
        g_atom = mol.GetAtomWithIdx(g_atom_idx)
        for p_atom_idx in phosphate_match:
            p_atom = mol.GetAtomWithIdx(p_atom_idx)
            for neighbor in g_atom.GetNeighbors():
                 if neighbor.GetIdx() == p_atom_idx:
                    glycerol_phosphate_atom = g_atom_idx
                    break
        if glycerol_phosphate_atom:
            break

    glycerol_ester_atom = None
    for e_match in ester_matches:
        ester_atom_idx = e_match[0]
        ester_atom = mol.GetAtomWithIdx(ester_atom_idx)
        for g_atom_idx in glycerol_match:
            g_atom = mol.GetAtomWithIdx(g_atom_idx)
            for neighbor in g_atom.GetNeighbors():
                if neighbor.GetIdx() == ester_atom_idx:
                    glycerol_ester_atom = g_atom_idx
                    break
            if glycerol_ester_atom:
                break
    
    if glycerol_phosphate_atom == glycerol_ester_atom:
        return False, "Phosphate and ester groups connected to the same carbon of the glycerol backbone."
     
    # Count rotatable bonds connected to the ester carbon and not the phosphate group - to ensure a single fatty acid chain.
    rotatable_bonds = 0
    for e_match in ester_matches:
        ester_c_atom_idx = e_match[1]
        ester_c_atom = mol.GetAtomWithIdx(ester_c_atom_idx)
        for neighbor in ester_c_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            is_phosphate_neighbor = False
            for p_atom_idx in phosphate_match:
                 if neighbor_idx == p_atom_idx:
                    is_phosphate_neighbor = True
                    break
            if not is_phosphate_neighbor:
                
               # Recursive search for rotatable bonds from the ester to count the fatty chain
                visited_atoms = set()
                stack = [neighbor_idx]
                while stack:
                    current_atom_idx = stack.pop()
                    if current_atom_idx in visited_atoms:
                        continue
                    visited_atoms.add(current_atom_idx)
                    current_atom = mol.GetAtomWithIdx(current_atom_idx)
                    
                    if current_atom.GetAtomicNum() == 6:
                        if current_atom.GetHybridization() != Chem.HybridizationType.SP:
                            rotatable_bonds += 1
                        for n in current_atom.GetNeighbors():
                            if n.GetIdx() not in visited_atoms:
                                stack.append(n.GetIdx())
    if rotatable_bonds < 2:
        return False, "Fatty acid chain is too short."
    
    # Check the number of phosphorus atoms - must be 1.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
         return False, f"Must have exactly 1 phosphorus atom, found {p_count}"

    
    return True, "Contains a glycerol backbone with a single fatty acid chain and a phosphate group"