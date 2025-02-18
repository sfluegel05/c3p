"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: neoflavonoid (CHEBI:XXXXX)
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is a 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define the coumarin (1-benzopyran-2-one) core pattern
    coumarin_core = Chem.MolFromSmarts('[O]=C1OC2=C(C=CC=C2)C=C1')
    core_matches = mol.GetSubstructMatches(coumarin_core)
    
    if not core_matches:
        return False, "No 1-benzopyran core found"
    
    # Check each core match for aryl substituent at position 4
    for match in core_matches:
        # Get the oxygen atom in the coumarin core
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:  # Oxygen atom
                # Find adjacent carbon in the pyrone ring (single bond)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.BondType.SINGLE for bond in atom.GetBondsToAtoms([neighbor.GetIdx()])):
                        # Check substituents on this carbon
                        for substituent in neighbor.GetNeighbors():
                            if substituent.GetIdx() in match:
                                continue  # Part of the core, skip
                            # Check if substituent is part of an aromatic ring
                            if substituent.GetIsAromatic():
                                return True, "1-benzopyran with aryl at position 4"
                            # Check if substituent is connected to an aromatic ring
                            elif substituent.GetDegree() > 0:
                                # Check neighboring atoms for aromatic rings
                                for sub_neighbor in substituent.GetNeighbors():
                                    if sub_neighbor.GetIsAromatic():
                                        return True, "1-benzopyran with aryl at position 4"
    
    return False, "No aryl substituent at position 4"