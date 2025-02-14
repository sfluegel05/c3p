"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:24842 indole alkaloid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is an alkaloid containing an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define indole SMARTS pattern
    indole_pattern = Chem.MolFromSmarts('c1c[cH]c2c1cccc2')  # Indole ring
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole skeleton found"

    # Identify basic nitrogen atoms (exclude aromatic nitrogens)
    basic_nitrogen_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            if not atom.GetIsAromatic():
                # Exclude amide nitrogens (N-C=O)
                is_amide = False
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:  # Carbon neighbor
                        for bond in mol.GetBonds():
                            if bond.GetBeginAtomIdx() == atom.GetIdx() and bond.GetEndAtomIdx() == neighbor.GetIdx():
                                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                    # Check if carbon is double-bonded to oxygen
                                    for nei_of_nei in neighbor.GetNeighbors():
                                        if nei_of_nei.GetAtomicNum() == 8:  # Oxygen
                                            for bond2 in mol.GetBonds():
                                                if bond2.GetBeginAtomIdx() == neighbor.GetIdx() and bond2.GetEndAtomIdx() == nei_of_nei.GetIdx():
                                                    if bond2.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                                        is_amide = True
                                                        break
                                        if is_amide:
                                            break
                            if is_amide:
                                break
                    if is_amide:
                        break
                if not is_amide:
                    basic_nitrogen_found = True
                    break  # Found at least one basic nitrogen atom

    if not basic_nitrogen_found:
        return False, "No basic nitrogen atoms found"

    return True, "Contains indole skeleton with at least one basic nitrogen atom (alkaloid)"