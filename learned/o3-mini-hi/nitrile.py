"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: A compound having the structure RC#N (nitrile)
Definition: A nitrile is a compound in which a carbon atom forms a triple bond with nitrogen,
with at least one carbon substituent attached to the nitrile carbon (i.e. it is a C-substituted derivative of HCN).
"""

from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is characterized by a carbon-nitrogen triple bond (C≡N) where the carbon is
    substituted with at least one non-hydrogen atom (i.e., RC#N rather than HC#N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if at least one nitrile (RC#N) group is found, False otherwise.
        str: A reason for the classification.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    nitrile_found = False
    reasons = []

    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        # Check if the bond is a triple bond
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Check for a bond between Carbon (atomic number 6) and Nitrogen (atomic number 7)
            # We consider only cases where carbon is the substituted atom.
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7:
                nitrile_carbon = atom1
                nitrile_nitrogen = atom2
            elif atom2.GetAtomicNum() == 6 and atom1.GetAtomicNum() == 7:
                nitrile_carbon = atom2
                nitrile_nitrogen = atom1
            else:
                continue  # Not a C≡N bond
            
            # Check that the nitrile carbon has at least one neighbor other than the nitrile nitrogen,
            # and that neighbor is not hydrogen (atomic number 1). This ensures it is a C-substituted nitrile.
            has_valid_substituent = False
            for neighbor in nitrile_carbon.GetNeighbors():
                if neighbor.GetIdx() != nitrile_nitrogen.GetIdx() and neighbor.GetAtomicNum() != 1:
                    has_valid_substituent = True
                    break
            
            if has_valid_substituent:
                nitrile_found = True
                reasons.append("Found C≡N group with substituted carbon")
    
    if nitrile_found:
        return True, "; ".join(reasons)
    
    return False, "No valid nitrile (RC#N) group found"