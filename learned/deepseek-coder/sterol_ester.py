"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:37664 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sterol backbone pattern (more general tetracyclic ring system)
    sterol_pattern = Chem.MolFromSmarts("[C@@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2C)CC=C4[C@@]3(CC[C@@H](C4)O)C")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"

    # Look for ester group (-O-C(=O)-) attached to the 3-hydroxy group
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester group found"

    # Find the oxygen atom of the 3-hydroxy group
    hydroxy_oxygen = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetTotalNumHs() == 1:
                    hydroxy_oxygen = atom
                    break
            if hydroxy_oxygen:
                break

    if not hydroxy_oxygen:
        return False, "No 3-hydroxy group found"

    # Check if any ester oxygen is attached to the 3-hydroxy carbon
    ester_attached = False
    for match in ester_matches:
        ester_oxygen = mol.GetAtomWithIdx(match[0])
        for neighbor in ester_oxygen.GetNeighbors():
            if neighbor.GetIdx() == hydroxy_oxygen.GetIdx():
                ester_attached = True
                break
        if ester_attached:
            break

    if not ester_attached:
        return False, "Ester group not attached to the 3-hydroxy group"

    return True, "Contains sterol backbone with ester linkage at the 3-hydroxy group"