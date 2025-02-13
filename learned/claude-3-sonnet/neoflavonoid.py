"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: CHEBI:75768 neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    Neoflavonoids have a 1-benzopyran core with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic benzopyran core pattern
    # More general pattern for the core structure
    benzopyran_pattern = Chem.MolFromSmarts("O1CCc2ccccc2C1")
    if benzopyran_pattern is None:
        return False, "Invalid benzopyran SMARTS pattern"
    
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No benzopyran core found"

    # Check for aromatic ring at position 4
    # This pattern looks for an aromatic ring connected at position 4
    aryl_pattern = Chem.MolFromSmarts("O1CCc2ccccc2C1c3ccccc3")
    if aryl_pattern is None:
        return False, "Invalid aryl SMARTS pattern"
    
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "No aryl substituent at position 4"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient number of rings"

    # Check for common substituents
    # Look for common oxygen-containing groups
    oxygen_pattern = Chem.MolFromSmarts("[OH,OC]")
    if oxygen_pattern is None:
        return False, "Invalid oxygen SMARTS pattern"
        
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    
    # Most neoflavonoids have oxygen-containing substituents
    if len(oxygen_matches) < 2:  # At least the pyran oxygen and one more
        return False, "Insufficient oxygen-containing substituents"

    # Verify aromatic character
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 6:
        return False, "Insufficient aromatic character"

    # Additional structural requirements
    # Check for presence of carbonyl group (often present in neoflavonoids)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    if carbonyl_pattern is not None and mol.HasSubstructMatch(carbonyl_pattern):
        confidence = "high"
    else:
        confidence = "medium"

    # Check molecular weight
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for neoflavonoid"

    # Count carbons to ensure reasonable size
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 15:
        return False, "Too few carbons for neoflavonoid structure"

    return True, f"Contains benzopyran core with aryl substituent at position 4 ({confidence} confidence)"