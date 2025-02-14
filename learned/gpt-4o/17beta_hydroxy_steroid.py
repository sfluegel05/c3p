"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid is defined by having a steroid backbone with a beta-configured hydroxy group at position 17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General steroid backbone pattern - focusing on four rings common in steroids
    steroid_pattern = Chem.MolFromSmarts("[#6]12CC[C@H]3[C@@H](CCC4=C3C=CC=C4)CC[C@@H]1C2")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No suitable steroid backbone found"

    # Look for the 17beta-hydroxy group configuration
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H]1CC[C@@H]2C[C@@H]3C[C@H](CC4=CC=CC=C34)[H]CC[C@@H]2C1")
    if mol.HasSubstructMatch(hydroxy_pattern):
        return True, "17beta-Hydroxy steroid structure identified"

    # Check for a general attachment of a beta-oriented hydroxy group at C17
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Find carbon atoms
            idx = atom.GetIdx()
            if any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors()):  # Check for connected Oxygen (hydroxy)
                # Verify stereochemistry implies beta orientation
                if atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    return True, "17beta-hydroxy group with identifiable stereochemistry found"

    return False, "No 17beta-hydroxy steroid configuration detected"