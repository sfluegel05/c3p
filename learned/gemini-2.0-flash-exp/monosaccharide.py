"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon atoms and no glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Less than 3 carbon atoms"

    # Check for aldehyde or ketone group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[H]")  # C(=O)H
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4]") # C(=O)C
    
    has_carbonyl = mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(ketone_pattern)
    if not has_carbonyl:
        return False, "No aldehyde or ketone group"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2][H]")  # -OH
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Count other oxygens not in a carbonyl or a hydroxyl group
    other_oxygens = 0
    for atom in mol.GetAtoms():
         if atom.GetAtomicNum() == 8:
              if not atom.HasSubstituent(Chem.AtomFromSmarts("[H]")) and not atom.IsInRing(): # Skip carbonyl oxygen and hydroxyl
                    other_oxygens +=1

    if other_oxygens > 1 : # Glycosidic bond should not be present
          return False, "Possible glycosidic linkage"
    

    # Check for at least 1 OH group for each C atom
    if hydroxyl_count < carbon_count - 1 :
        return False, "Not enough hydroxyl groups"
    
    # Check for molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for monosaccharide"

    
    return True, "Meets criteria for a monosaccharide"