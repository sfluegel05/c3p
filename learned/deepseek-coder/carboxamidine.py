"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:37671 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine has the structure RC(=NR)NR2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxamidine pattern: RC(=NR)NR2
    # More specific pattern to match exactly one central carbon with =NR and NR2 groups
    carboxamidine_pattern = Chem.MolFromSmarts("[CX3](=[NX2])[NX3H2,NX3H1,NX3H0]")
    
    # Find all matches of the carboxamidine pattern
    matches = mol.GetSubstructMatches(carboxamidine_pattern)
    if not matches:
        return False, "No carboxamidine structure (RC(=NR)NR2) found"

    # Check each potential carboxamidine center
    for match in matches:
        central_carbon = match[0]
        double_bond_nitrogen = match[1]
        single_bond_nitrogen = match[2]
        
        # Verify the central carbon has exactly 3 bonds
        if len(mol.GetAtomWithIdx(central_carbon).GetBonds()) != 3:
            continue
            
        # Verify the double bond nitrogen has exactly 2 bonds
        if len(mol.GetAtomWithIdx(double_bond_nitrogen).GetBonds()) != 2:
            continue
            
        # Verify the single bond nitrogen has exactly 3 bonds
        if len(mol.GetAtomWithIdx(single_bond_nitrogen).GetBonds()) != 3:
            continue
            
        # Check that the single bond nitrogen is not part of a guanidine group
        # by ensuring it's not bonded to another carbon with =NR group
        guanidine_check = False
        for neighbor in mol.GetAtomWithIdx(single_bond_nitrogen).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for bond in neighbor.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        if bond.GetBeginAtom().GetAtomicNum() == 7 or bond.GetEndAtom().GetAtomicNum() == 7:
                            guanidine_check = True
                            break
                if guanidine_check:
                    break
                    
        if not guanidine_check:
            return True, "Contains the carboxamidine structure (RC(=NR)NR2)"

    return False, "No valid carboxamidine structure found after detailed checks"