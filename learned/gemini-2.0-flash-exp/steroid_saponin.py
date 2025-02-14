"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is a saponin derived from a hydroxysteroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for steroid core and glycosidic linkages
    # Improved steroid core pattern: four fused rings with any atoms and more flexibility
    steroid_core_pattern = Chem.MolFromSmarts("[*]12[*]34[*]([*](1)[*]2)[*]3[*](4)")
    if steroid_core_pattern is None:
        return None, "Could not create steroid core pattern" # sanity check
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found."
    
    # Check for at least one alcohol group attached to the core 
    hydroxy_group_pattern = Chem.MolFromSmarts("[*]~[OH]")
    if hydroxy_group_pattern is None:
        return None, "Could not create alcohol pattern" # sanity check
    if not mol.HasSubstructMatch(hydroxy_group_pattern):
         return False, "No hydroxy group found. Thus not a steroid saponin"
         
    # Check for glycosidic linkages (ring carbon-O-sugar carbon)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[C;R][O][C;R](O)")
    if glycosidic_linkage_pattern is None:
         return None, "Could not create glycosidic linkage pattern"  # sanity check
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    if len(glycosidic_matches) == 0 :
      return False, "No glycosidic linkage found. Thus not a saponin"
    
    # Check for sugar presence (ring containing multiple O)
    sugar_pattern = Chem.MolFromSmarts("[C;R]1[C;R](O)[C;R](O)[C;R]*[C;R]*[O;R]1")
    if sugar_pattern is None:
        return None, "Could not create sugar pattern" # sanity check
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0 :
      return False, "No sugar moiety found. Thus not a saponin"
    
    return True, "Contains steroid core with glycosidic linkages and sugar moieties (i.e. is a steroid saponin)."