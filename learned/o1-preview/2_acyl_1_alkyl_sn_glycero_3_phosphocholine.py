"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    A 2-acyl-1-alkyl-sn-glycero-3-phosphocholine has a glycerol backbone with an ether-linked alkyl chain at position 1,
    an ester-linked acyl chain at position 2, and a phosphocholine group at position 3.
    
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

    # Add hydrogens to preserve stereochemistry
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns for the required substructures

    # Glycerol backbone with positions labeled
    glycerol_pattern = Chem.MolFromSmarts("[C:1][C@H:2][C:3]")

    # Ether-linked alkyl chain at position 1 (carbon connected to oxygen connected to alkyl chain)
    ether_alkyl_pattern = Chem.MolFromSmarts("[C:1]-O-[C;R0]")

    # Ester-linked acyl chain at position 2
    ester_acyl_pattern = Chem.MolFromSmarts("[C@H:2]-O-C(=O)-[C;R0]")

    # Phosphocholine group at position 3
    phosphocholine_pattern = Chem.MolFromSmarts("[C:3]-O-P(=O)([O-])OCC[N+](C)(C)C")

    # Search for glycerol backbone
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone with correct stereochemistry found"

    # Iterate over matches to check substitutions
    for match in matches:
        c1_idx, c2_idx, c3_idx = match

        # Check for ether-linked alkyl chain at position 1
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        is_ether_alkyl = False
        for neighbor in c1_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for n_neighbor in neighbor.GetNeighbors():
                    if n_neighbor.GetIdx() != c1_idx and n_neighbor.GetAtomicNum() == 6 and not n_neighbor.IsInRing():
                        is_ether_alkyl = True
                        break
                if is_ether_alkyl:
                    break
        if not is_ether_alkyl:
            continue

        # Check for ester-linked acyl chain at position 2
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        is_ester_acyl = False
        for neighbor in c2_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for n_neighbor in neighbor.GetNeighbors():
                    if n_neighbor.GetIdx() != c2_idx and n_neighbor.GetAtomicNum() == 6:
                        # Check for carbonyl group
                        c_atom = n_neighbor
                        has_carbonyl = False
                        for c_neighbor in c_atom.GetNeighbors():
                            if c_neighbor.GetAtomicNum() == 8 and c_atom.GetBondBetweenAtoms(c_atom.GetIdx(), c_neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                has_carbonyl = True
                                break
                        if has_carbonyl:
                            # Check for attached alkyl chain
                            for c_neighbor in c_atom.GetNeighbors():
                                if c_neighbor.GetAtomicNum() == 6 and c_neighbor.GetIdx() != neighbor.GetIdx():
                                    is_ester_acyl = True
                                    break
                    if is_ester_acyl:
                        break
            if is_ester_acyl:
                break
        if not is_ester_acyl:
            continue

        # Check for phosphocholine group at position 3
        c3_atom = mol.GetAtomWithIdx(c3_idx)
        is_phosphocholine = False
        for neighbor in c3_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for n_neighbor in neighbor.GetNeighbors():
                    if n_neighbor.GetAtomicNum() == 15:  # Phosphorus
                        # Check for choline part
                        p_atom = n_neighbor
                        for p_neighbor in p_atom.GetNeighbors():
                            if p_neighbor.GetAtomicNum() == 8 and p_neighbor.GetIdx() != neighbor.GetIdx():
                                for o_neighbor in p_neighbor.GetNeighbors():
                                    if o_neighbor.GetAtomicNum() == 6:
                                        # Check for choline chain OCC[N+](C)(C)C
                                        choline_pattern = Chem.MolFromSmarts("OCC[N+](C)(C)C")
                                        if mol.HasSubstructMatch(choline_pattern):
                                            is_phosphocholine = True
                                            break
                            if is_phosphocholine:
                                break
                    if is_phosphocholine:
                        break
            if is_phosphocholine:
                break
        if not is_phosphocholine:
            continue

        # Check stereochemistry at C2 (should be specified)
        # Get chiral tag of C2 atom
        chiral_tag = c2_atom.GetChiralTag()
        if chiral_tag == Chem.ChiralType.CHI_UNSPECIFIED:
            return False, "C2 atom does not have specified stereochemistry"
        
        # If all checks pass
        return True, "Molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"

    return False, "Molecule does not match the required structure"