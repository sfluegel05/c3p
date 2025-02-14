"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define SMARTS patterns for common steroid cores with some flexibility
    steroid_core_patterns = [
        Chem.MolFromSmarts("[*]12[*]3[*]4[*]1[*]C[*]2[*]34"), # basic four fused rings
        Chem.MolFromSmarts("[*]12[*]3[*]4[*](=[*])1[*]C[*]2[*]34"), # with a double bond
        Chem.MolFromSmarts("[*]12[*]3[*]4[*]1[*]C[*](=[*])2[*]34"), # with a double bond at a different position
        Chem.MolFromSmarts("[*]12[*]3[*]4[*]1[*]C[*]2[*](=[*])34"),  # with a double bond at a different position
        Chem.MolFromSmarts("[*]12[*]3[*]4[*]1[*](=[*])C[*]2[*]34"),
        Chem.MolFromSmarts("[*]12[C]([C])([*])[*]3[*]4[*]1[*]C[*]2[*]34"), # with a methyl group
        Chem.MolFromSmarts("[*]12[C]([C])([*])[*]3[*]4[*](=[*])1[*]C[*]2[*]34"),
        Chem.MolFromSmarts("[*]12[C]([C])([*])[*]3[*]4[*](=[*])1[*]C[*](=[*])2[*]34"),
        Chem.MolFromSmarts("[*]12[C]([C])([*])[*]3[*]4[*](=[*])1[*](=[*])C[*]2[*]34"),
    ]
    
    has_steroid_core = False
    for pattern in steroid_core_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_steroid_core = True
            break
    if not has_steroid_core:
       return False, "No steroid core found"


    # Check for at least one alcohol group attached to the steroid core
    # Find steroid core matches
    core_matches = []
    for pattern in steroid_core_patterns:
      if pattern is not None:
        core_matches.extend(mol.GetSubstructMatches(pattern))

    has_hydroxy_on_core = False
    
    for match in core_matches:
        for atom_idx in match:
            # Check if the atoms in the core have an O attached by checking the neighbors of the atom indices
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
               if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1: # Check if neighbor is oxygen and has one H
                   has_hydroxy_on_core = True
                   break
            if has_hydroxy_on_core:
              break # move to next steroid match if one found
        if has_hydroxy_on_core:
            break # exit the loop after finding one hydroxyl group
            
    if not has_hydroxy_on_core:
       return False, "No hydroxyl group found on steroid core"
    
    # Check for glycosidic linkages (ring carbon-O-sugar carbon)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[C;R][O][C;R](O)")
    if glycosidic_linkage_pattern is None:
         return None, "Could not create glycosidic linkage pattern"
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    if len(glycosidic_matches) == 0 :
      return False, "No glycosidic linkage found. Thus not a saponin"

    # Check for sugar presence (ring containing multiple O)
    sugar_pattern = Chem.MolFromSmarts("[C;R]1[C;R](O)[C;R](O)[C;R]*[C;R]*[O;R]1")
    if sugar_pattern is None:
        return None, "Could not create sugar pattern"
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0 :
      return False, "No sugar moiety found. Thus not a saponin"


    return True, "Contains steroid core with glycosidic linkages and sugar moieties (i.e. is a steroid saponin)."