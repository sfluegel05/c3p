"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is characterized by a steroid structure with a hydroxyl
    group at the 3rd carbon in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone: expanded pattern to consider multiple conditions
    # Recognizing a tetracyclic ABCD ring framework common in steroids
    steroid_pattern = Chem.MolFromSmarts("C1([C@H]2[C@H]([C@H]3CC[C@@H]4C=C/[C@@H]([C@]4([H])[H])CCCC3)CCC2)C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define 3beta-hydroxy group: Consider more generic position for hydroxyl
    hydroxy_attributes = {
        'Carbon': (6, ),  # Only Carbon atoms
        'Hydroxy': Chem.MolFromSmarts("[OH]"),  # Hydroxy group
    }

    # Find 3beta-hydroxy group
    beta_hydroxy = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6: # carbon
            # Check if adjacent atom is an OH
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8: # oxygen
                    beta_hydroxy = True
                    break
        if beta_hydroxy:
            break

    if not beta_hydroxy:
        return False, "No 3beta-hydroxy group found"
    
    return True, "Contains steroid backbone with 3beta-hydroxy group"