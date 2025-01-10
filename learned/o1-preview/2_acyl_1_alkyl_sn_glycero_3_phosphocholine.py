"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: CHEBI:63866 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    This class has:
    - A glycerol backbone.
    - An ether-linked alkyl chain at position 1.
    - An ester-linked acyl chain at position 2.
    - A phosphocholine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to ensure valence consistency
    try:
        Chem.Kekulize(mol)
    except:
        pass

    # Define SMARTS patterns for the required substructures

    # Glycerol backbone without stereochemistry specification
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Ether-linked alkyl chain at position 1
    # An ether linkage (C-O-C) where one of the carbons is part of the glycerol backbone
    ether_pattern = Chem.MolFromSmarts("COCC(O)CO")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No ether-linked alkyl chain at position 1 found"

    # Ester-linked acyl chain at position 2
    # An ester linkage (O-C=O) where the oxygen is connected to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ester_connected = False
    for match in ester_matches:
        ester_oxygen_idx = match[0]
        glycerol_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)
        for bond in glycerol_oxygen.GetBonds():
            neighbor = bond.GetOtherAtom(glycerol_oxygen)
            # Check if the oxygen is connected to a carbon of the glycerol backbone
            if neighbor.GetAtomicNum() == 6:
                paths = Chem.rdmolops.GetShortestPath(mol, neighbor.GetIdx(), ester_oxygen_idx)
                if len(paths) <= 3:
                    ester_connected = True
                    break
        if ester_connected:
            break
    if not ester_connected:
        return False, "No ester-linked acyl chain at position 2 found"

    # Phosphocholine group at position 3
    # Phosphate group with a choline moiety
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)(OCC[N+](C)(C)C)[O-]")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group at position 3 found"

    # Additional checks to ensure correct connectivity
    # Find glycerol backbone carbons
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    for match in glycerol_matches:
        # Indices of the glycerol carbons and oxygens
        o1_idx = match[0]
        c1_idx = match[1]
        o2_idx = match[2]
        c2_idx = match[3]
        o3_idx = match[4]

        # Check for ether linkage at position 1 (connected to c1)
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        ether_found = False
        for bond in c1_atom.GetBonds():
            nbr = bond.GetOtherAtom(c1_atom)
            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != o1_idx:
                # Oxygen connected to an alkyl chain
                for obond in nbr.GetBonds():
                    onbr = obond.GetOtherAtom(nbr)
                    if onbr.GetAtomicNum() == 6 and onbr.GetIdx() != c1_idx:
                        ether_found = True
                        break
            if ether_found:
                break
        if not ether_found:
            continue  # Try next glycerol match

        # Check for ester linkage at position 2 (connected to o2)
        o2_atom = mol.GetAtomWithIdx(o2_idx)
        ester_found = False
        for bond in o2_atom.GetBonds():
            nbr = bond.GetOtherAtom(o2_atom)
            if nbr.GetAtomicNum() == 6:
                # Carbon connected to carbonyl oxygen
                carbonyl_found = False
                for cbond in nbr.GetBonds():
                    cnbr = cbond.GetOtherAtom(nbr)
                    if cnbr.GetAtomicNum() == 8 and cbond.GetBondTypeAsDouble() == 2.0:
                        carbonyl_found = True
                        break
                if carbonyl_found:
                    ester_found = True
                    break
        if not ester_found:
            continue  # Try next glycerol match

        # Check for phosphocholine group at position 3 (connected to o3)
        o3_atom = mol.GetAtomWithIdx(o3_idx)
        phosphate_found = False
        for bond in o3_atom.GetBonds():
            nbr = bond.GetOtherAtom(o3_atom)
            if nbr.GetAtomicNum() == 15:  # Phosphorus
                # Check for choline moiety
                for pbond in nbr.GetBonds():
                    pnbr = pbond.GetOtherAtom(nbr)
                    if pnbr.GetAtomicNum() == 8 and pbond.GetBondTypeAsDouble() == 1.0:
                        # Oxygen connected to choline
                        for obond in pnbr.GetBonds():
                            onbr = obond.GetOtherAtom(pnbr)
                            if onbr.GetAtomicNum() == 6:
                                # Carbon chain leading to nitrogen
                                path = Chem.rdmolops.GetShortestPath(mol, onbr.GetIdx(), o3_idx)
                                atoms_in_path = [mol.GetAtomWithIdx(idx) for idx in path]
                                nitrogen_found = any(atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1 for atom in atoms_in_path)
                                if nitrogen_found:
                                    phosphate_found = True
                                    break
                    if phosphate_found:
                        break
            if phosphate_found:
                break
        if not phosphate_found:
            continue  # Try next glycerol match

        # All substructures found with correct connectivity
        return True, "Molecule matches all required substructures for 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"

    return False, "Molecule does not match all required substructures"