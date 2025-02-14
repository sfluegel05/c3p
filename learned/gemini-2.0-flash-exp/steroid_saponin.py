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
    # Simplified steroid core pattern: four fused rings with some variability.
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C][C]3[C][C]4[C]([C]1)[C]2[C]34")
    if steroid_core_pattern is None:
        return None, "Could not create steroid core pattern" # sanity check
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found."
    
    # Check for glycosidic linkages (ring carbon-O-sugar carbon)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[C;R][O][C;R](O)")
    if glycosidic_linkage_pattern is None:
         return None, "Could not create glycosidic linkage pattern"  # sanity check
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    if len(glycosidic_matches) == 0 :
      return False, "No glycosidic linkage found. Thus not a saponin"
    
    # Check for sugar presence (ring containing multiple O)
    sugar_pattern = Chem.MolFromSmarts("[C]1[C](O)[C](O)[C](O)[C](O)[C]1")
    if sugar_pattern is None:
        return None, "Could not create sugar pattern" # sanity check
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0 :
      return False, "No sugar moiety found. Thus not a saponin"
    
    return True, "Contains steroid core with glycosidic linkages and sugar moieties (i.e. is a steroid saponin)."