"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for basic steroid core (four rings)
    steroid_core = Chem.MolFromSmarts("C1CC2CCC3C4CCCC4CCC3C2C1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"
    
    # Look for carboxylic acid group or its conjugates
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    glycine_conjugate = Chem.MolFromSmarts("[CX3](=O)NCC(=O)[OH]")
    taurine_conjugate = Chem.MolFromSmarts("[CX3](=O)NCCS(=O)(=O)[OH]")
    
    has_acid = mol.HasSubstructMatch(carboxylic_acid)
    has_glycine = mol.HasSubstructMatch(glycine_conjugate)
    has_taurine = mol.HasSubstructMatch(taurine_conjugate)
    
    if not (has_acid or has_glycine or has_taurine):
        return False, "No carboxylic acid group or conjugates found"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 27:  # Most bile acids have 24 carbons
        return False, f"Carbon count ({c_count}) outside typical range for bile acids"
    
    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyls == 0:
        return False, "No hydroxyl groups found"
        
    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 700:
        return False, f"Molecular weight ({mol_wt}) outside typical range for bile acids"
    
    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"
    
    # Check for reasonable number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too few rotatable bonds"
        
    # Check for sp3 carbons (should have many in steroid core)
    sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if sp3_carbons < 15:
        return False, "Too few sp3 carbons for bile acid structure"

    return True, "Matches bile acid structure with steroid core, appropriate functional groups, and typical molecular features"