"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flag to indicate if phospho sugar is found
    phospho_sugar_found = False

    # Iterate over all phosphorus atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Atomic number of phosphorus
            p_atom = atom
            neighbors = p_atom.GetNeighbors()
            double_bonded_oxygens = 0
            single_bonded_oxygens = []
            # Analyze the bonds to phosphorus
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_bonded_oxygens += 1
                    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        single_bonded_oxygens.append(neighbor)
            # Check for phosphate group: P bonded to one double-bonded O and three single-bonded O
            if double_bonded_oxygens == 1 and len(single_bonded_oxygens) >=2:
                # Check for ester linkage to carbon
                for o_atom in single_bonded_oxygens:
                    # Check if oxygen is connected to carbon
                    for neighbor in o_atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != p_atom.GetIdx():
                            c_atom = neighbor
                            # Now check if carbon is part of a sugar ring
                            # Get rings that include this carbon
                            ring_info = mol.GetRingInfo()
                            atom_rings = ring_info.AtomRings()
                            for ring in atom_rings:
                                if len(ring) == 5 or len(ring) == 6:
                                    if c_atom.GetIdx() in ring:
                                        # Check if ring consists of C and O atoms
                                        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                                        ring_atom_nums = [atom.GetAtomicNum() for atom in ring_atoms]
                                        if all(num == 6 or num == 8 for num in ring_atom_nums):
                                            # Check if carbons in ring have hydroxyl groups (O-H)
                                            has_hydroxyls = True
                                            for ring_atom in ring_atoms:
                                                if ring_atom.GetAtomicNum() == 6:
                                                    hydroxyl_found = False
                                                    for neighbor in ring_atom.GetNeighbors():
                                                        if neighbor.GetAtomicNum() == 8:
                                                            bond = mol.GetBondBetweenAtoms(ring_atom.GetIdx(), neighbor.GetIdx())
                                                            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                                                # Check if oxygen has a hydrogen (is a hydroxyl group)
                                                                if neighbor.GetTotalNumHs() > 0:
                                                                    hydroxyl_found = True
                                                                    break
                                                    if not hydroxyl_found:
                                                        has_hydroxyls = False
                                                        break
                                            if has_hydroxyls:
                                                phospho_sugar_found = True
                                                return True, "Contains monosaccharide with an alcoholic hydroxy group esterified with phosphoric acid"
    if not phospho_sugar_found:
        return False, "No phospho sugar moiety found"