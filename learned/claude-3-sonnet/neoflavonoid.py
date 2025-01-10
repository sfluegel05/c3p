"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: CHEBI:75768 neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

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

    # Basic benzopyran core pattern - allowing for various bond types
    # [#6] represents any carbon, ~ represents any bond
    benzopyran_pattern = Chem.MolFromSmarts("O1[#6]~[#6]c2ccccc2[#6]1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No benzopyran core found"

    # Check for aromatic ring at position 4
    # More flexible pattern for aryl substituent
    aryl_pattern = Chem.MolFromSmarts("O1[#6]~[#6]c2ccccc2[#6]1[#6]3[#6]~[#6]~[#6]~[#6]~[#6]3")
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "No aryl substituent at position 4"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient number of rings"

    # Check for common substituents
    # Look for oxygen-containing groups (hydroxyl, methoxy, carbonyl)
    oxygen_pattern = Chem.MolFromSmarts("[$([OH]),$([O][CH3]),$(O=[#6])]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    
    # Most neoflavonoids have multiple oxygen-containing substituents
    if len(oxygen_matches) < 1:  # At least one besides the pyran oxygen
        return False, "Insufficient oxygen-containing substituents"

    # Verify aromatic character
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 6:
        return False, "Insufficient aromatic character"

    # Check molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside typical range for neoflavonoids"

    # Count carbons to ensure reasonable size
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 15:
        return False, "Too few carbons for neoflavonoid structure"

    # Additional checks for common neoflavonoid features
    confidence = "medium"
    
    # Check for common substitution patterns
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    methoxy_pattern = Chem.MolFromSmarts("[O][CH3]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    
    if (mol.HasSubstructMatch(hydroxy_pattern) and 
        mol.HasSubstructMatch(carbonyl_pattern)):
        confidence = "high"
    
    # Check for common positions of substitution
    if mol.HasSubstructMatch(Chem.MolFromSmarts("O1[#6]~[#6]c2c(O)cccc2[#6]1")):
        confidence = "high"

    return True, f"Contains benzopyran core with aryl substituent at position 4 ({confidence} confidence)"