"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a molecule containing a benzene ring with four chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "Molecule does not contain a benzene ring"
    
    # 2. Check for exactly four chlorine atoms
    chlorine_count = 0
    chlorine_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:
            chlorine_count += 1
            chlorine_atoms.append(atom)
    if chlorine_count != 4:
        return False, f"Molecule does not have exactly 4 chlorine atoms, found {chlorine_count}"
    
    # 3. Check if all four chlorines are linked (directly or via one atom) to benzene ring
    for chlorine in chlorine_atoms:
        is_linked_to_benzene = False
        for neighbor in chlorine.GetNeighbors():
            if neighbor.HasSubstructMatch(Chem.MolFromSmarts("[c]")):
                 is_linked_to_benzene = True
                 break
            else:
                for neighbor2 in neighbor.GetNeighbors():
                    if neighbor2.HasSubstructMatch(Chem.MolFromSmarts("[c]")):
                       is_linked_to_benzene = True
                       break
                if is_linked_to_benzene:
                  break
        if not is_linked_to_benzene:
            return False, "Not all chlorine atoms are directly or indirectly linked to benzene ring"


    return True, "Molecule contains a benzene ring with exactly four chlorine substituents"