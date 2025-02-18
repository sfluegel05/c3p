"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: CHEBI:27306 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D analog based on its SMILES string.
    Vitamin D compounds are seco-steroids with a broken B ring and conjugated triene system.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a vitamin D analog, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for at least one hydroxyl group
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 1:
        return False, "No hydroxyl groups"
    
    # Check for conjugated triene system (three conjugated double bonds)
    triene_pattern = Chem.MolFromSmarts('[CH1]=[CH1]-[CH1]=[CH1]-[CH1]=[CH1]')
    if not mol.HasSubstructMatch(triene_pattern):
        # Alternative pattern for conjugated triene with possible branching
        triene_pattern_flex = Chem.MolFromSmarts('*~*=*~*=*~*')
        if not mol.HasSubstructMatch(triene_pattern_flex):
            return False, "No conjugated triene system"
    
    # Check for seco-steroid core (simplified as presence of four rings, one of which is open)
    # This is a heuristic and may not cover all cases
    ri = mol.GetRingInfo()
    if len(ri.AtomRings()) < 3:  # Original steroids have 4 rings, seco may have 3
        return False, "Insufficient rings for seco-steroid structure"
    
    # Check molecular weight to filter out small molecules (vitamin D analogs are typically >350 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 350:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"
    
    return True, "Has hydroxyl groups, conjugated triene, and seco-steroid features"