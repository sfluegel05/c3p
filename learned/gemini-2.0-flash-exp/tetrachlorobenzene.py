"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a molecule containing a benzene ring with four chlorine substituents (directly or via one atom).

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

    # Check for a benzene ring with four chlorine substituents
    # This SMARTS pattern looks for a benzene ring with at least four chlorine atoms attached directly or via one atom,
    # the X means any atom.
    
    tetrachloro_benzene_pattern = Chem.MolFromSmarts("c1([X,Cl])([X,Cl])([X,Cl])([X,Cl])cc1")
    if not mol.HasSubstructMatch(tetrachloro_benzene_pattern):
      
        # If the above fails check a similar pattern in which the chlorines are attached to the ring
        tetrachloro_benzene_pattern_direct = Chem.MolFromSmarts("c1(Cl)(Cl)(Cl)(Cl)cc1")
        if not mol.HasSubstructMatch(tetrachloro_benzene_pattern_direct):
          
            # If the above fails, return that it does not contain a tetrachlorobenzene
            return False, "Molecule does not contain a benzene ring with four chlorine substituents"

    #Count the atoms to make sure only 4 are found in total, this is an extra safety check
    chlorine_count = 0
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() == 17:
        if atom.HasSubstructMatch(Chem.MolFromSmarts("[c]")):
          chlorine_count += 1
        else:
          for neighbor in atom.GetNeighbors():
            if neighbor.HasSubstructMatch(Chem.MolFromSmarts("[c]")):
              chlorine_count += 1


    if chlorine_count != 4:
        return False, f"Found {chlorine_count} chlorine atoms on the benzene ring, need exactly 4"

    return True, "Molecule contains a benzene ring with exactly four chlorine substituents"