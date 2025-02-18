"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI: beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate has a beta-D-glucuronic acid moiety (deprotonated carboxylate)
    linked via an O-glycosidic bond with correct stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define beta-D-glucuronide pattern with atom map on anomeric carbon
    # Pattern: carboxylate connected to C4 (D-configuration), beta linkage at C1
    # SMARTS accounts for beta configuration (axial C1-O-C-glycone)
    glucuronide_pattern = Chem.MolFromSmarts(
        "[O-]C(=O)[C@@H:1]1O[C@H]([C@H]([C@@H]([C@H]1O)O)O)O"
    )
    if not glucuronide_pattern:
        return False, "Invalid pattern"
    
    matches = mol.GetSubstructMatches(glucuronide_pattern)
    if not matches:
        return False, "No beta-D-glucuronide substructure"
    
    # Check each match for glycosidic bond from anomeric carbon
    for match in matches:
        # Find the anomeric carbon (atom map 1 in pattern)
        anomeric_idx = None
        for i, atom in enumerate(glucuronide_pattern.GetAtoms()):
            if atom.GetAtomMapNum() == 1:
                anomeric_idx = i
                break
        if anomeric_idx is None:
            continue  # Shouldn't happen if pattern is correct
        
        anomeric_carbon = mol.GetAtomWithIdx(match[anomeric_idx])
        
        # Check for oxygen neighbor connected to non-glucuronide atoms
        for neighbor in anomeric_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                # Check if oxygen is bonded to atom outside the glucuronide match
                for bond in neighbor.GetBonds():
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetIdx() not in match:
                        return True, "Contains beta-D-glucuronide with O-glycosidic bond"
    
    return False, "No valid O-glycosidic bond in glucuronide structure"