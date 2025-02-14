"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene has a benzene ring with exactly three chlorine atoms attached to it
    either directly or via a single bond to another atom (i.e. -O-Cl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total number of chlorine atoms
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    if len(chlorine_atoms) != 3:
        return False, f"Found {len(chlorine_atoms)} chlorine atoms, need exactly 3"


    # Check for benzene ring with attached chlorines, either directly or via single bond.
    # "[cX3]([Cl])" matches a benzene carbon connected to a Cl atom.
    # "[cX3]-[#6,#7,#8,#15,#16]-[Cl]" matches a benzene carbon connected to an atom
    # that is not hydrogen, then connected to a Cl atom.
    pattern1 = Chem.MolFromSmarts("c1([cX3]([Cl]))ccccc1")
    pattern2 = Chem.MolFromSmarts("c1([cX3]-[#6,#7,#8,#15,#16]-[Cl])ccccc1")
    
    if mol.HasSubstructMatch(pattern1):
          count = len(mol.GetSubstructMatches(pattern1))

          
          if count == 1 and not mol.HasSubstructMatch(pattern2):
            return True, "Contains a benzene ring with exactly 3 chlorine atoms attached"
          elif count == 1 and mol.HasSubstructMatch(pattern2):
                pattern3 = Chem.MolFromSmarts("c1([cX3]([Cl]))([cX3]-[#6,#7,#8,#15,#16]-[Cl])ccccc1")
                pattern4 = Chem.MolFromSmarts("c1([cX3]-[#6,#7,#8,#15,#16]-[Cl])([cX3]-[#6,#7,#8,#15,#16]-[Cl])ccccc1")
                if not (mol.HasSubstructMatch(pattern3) or mol.HasSubstructMatch(pattern4)):
                     return True, "Contains a benzene ring with exactly 3 chlorine atoms attached"
          
          elif count == 2:
               pattern3 = Chem.MolFromSmarts("c1([cX3]([Cl]))([cX3]([Cl]))ccccc1")
               if mol.HasSubstructMatch(pattern3):
                    pattern5 = Chem.MolFromSmarts("c1([cX3]([Cl]))([cX3]([Cl]))([cX3]-[#6,#7,#8,#15,#16]-[Cl])ccccc1")
                    if not mol.HasSubstructMatch(pattern5):
                        return True, "Contains a benzene ring with exactly 3 chlorine atoms attached"
          elif count == 3:
                pattern5 = Chem.MolFromSmarts("c1([cX3]([Cl]))([cX3]([Cl]))([cX3]([Cl]))ccccc1")
                if mol.HasSubstructMatch(pattern5):
                    return True, "Contains a benzene ring with exactly 3 chlorine atoms attached"
               


    if mol.HasSubstructMatch(pattern2):
        count2 = len(mol.GetSubstructMatches(pattern2))
        if count2==1:
             return False, "Found 1 chlorine atoms attached to benzene via a single bond, need exactly 3"
        elif count2==2:
              pattern6 = Chem.MolFromSmarts("c1([cX3]-[#6,#7,#8,#15,#16]-[Cl])([cX3]-[#6,#7,#8,#15,#16]-[Cl])ccccc1")
              if mol.HasSubstructMatch(pattern6):
                pattern7 = Chem.MolFromSmarts("c1([cX3]-[#6,#7,#8,#15,#16]-[Cl])([cX3]-[#6,#7,#8,#15,#16]-[Cl])([cX3]([Cl]))ccccc1")
                if not mol.HasSubstructMatch(pattern7):
                     return True, "Contains a benzene ring with exactly 3 chlorine atoms attached"
        elif count2 == 3:
              pattern8 = Chem.MolFromSmarts("c1([cX3]-[#6,#7,#8,#15,#16]-[Cl])([cX3]-[#6,#7,#8,#15,#16]-[Cl])([cX3]-[#6,#7,#8,#15,#16]-[Cl])ccccc1")
              if mol.HasSubstructMatch(pattern8):
                   return True, "Contains a benzene ring with exactly 3 chlorine atoms attached"


    
    return False, "Did not find a benzene with exactly 3 chlorine atoms attached directly or via single bond."