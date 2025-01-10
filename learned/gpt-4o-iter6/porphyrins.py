"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins have a macrocyclic structure consisting of four pyrrole rings connected by methine bridges,
    potentially coordinated with metal ions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the core porphyrin macrocycle
    # Porphyrins: [N]!@[C]1=[C]ncc1-[C]2=[C]ncc2-[C]3=[C]ncc3-[C]4=[C]ncc4
    # Flexibility added with variation in pyrrole ring possible links
    porphyrin_pattern = Chem.MolFromSmarts("[nR1]1[cR2][cR2][cR2]2[nR1][cR2][cR2]3[cR2][nR1][cR2][cR2]4[cR2][nR1][cR2][cR2]1[cR2]4[cR2]32")

    if porphyrin_pattern is None:
        return None, "Error in SMARTS pattern definition"

    # Attempt to find the macrocyclic structure
    if mol.HasSubstructMatch(porphyrin_pattern):
        # Note: checking metal coordination could be checked by looking for metal atoms in connectivity
        metal_involved = any(atom.GetAtomicNum() in [24, 25, 26, 27, 28, 29, 30, 12] for atom in mol.GetAtoms())  # common metals like Mg, Fe, etc.
        return True, f"SMILES string contains a porphyrin-like macrocyclic structure, {'with metal coordination' if metal_involved else 'without metal coordination'}"
    else:
        return False, "No macrocyclic structure typical of porphyrins found in SMILES string"