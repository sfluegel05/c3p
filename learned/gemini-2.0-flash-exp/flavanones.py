"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones have a 2-phenyl-3,4-dihydro-2H-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for the flavanone core structure
    # Reduced the specificity of SMARTS to allow for substitution on the rings
    flavanone_core_smarts = "[cH]1[cH][cH][c]2[O][CH]([CH2][CH]2[C]=O)[c]1"
    flavanone_core_pattern = Chem.MolFromSmarts(flavanone_core_smarts)
    
    if not mol.HasSubstructMatch(flavanone_core_pattern):
       return False, "Flavanone core structure not found"
    
    
    # check for a phenyl group at the 2 position
    phenyl_smarts = "[cH]1[cH][cH][cH][cH][cH]1"
    phenyl_pattern = Chem.MolFromSmarts(phenyl_smarts)
    
    matches = mol.GetSubstructMatches(flavanone_core_pattern)
    
    if matches:
      match = matches[0]
      for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 1: # This corresponds to the carbon on the 2-position of the chromane core

                  for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() == 0:  # This corresponds to the carbon from the C in phenyl

                        attached_phenyl = False
                        for phenyl_neighbor in neighbor.GetNeighbors():
                          if phenyl_neighbor.GetSymbol() == 'C' and phenyl_neighbor.IsInRing() and phenyl_neighbor.GetTotalNumHs() == 1:
                              attached_phenyl = True
                              break
                        if not attached_phenyl:
                            return False, "Phenyl substituent not found at correct position"

                        
                        #Check for the phenyl group using substruct match on the whole molecule
                        if not mol.HasSubstructMatch(phenyl_pattern):
                          return False, "Phenyl substituent not found"

                        break
                  break
    
    
    return True, "Molecule is a flavanone"