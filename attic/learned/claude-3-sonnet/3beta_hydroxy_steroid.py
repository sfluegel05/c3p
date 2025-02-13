"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for basic steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for 3-OH group with specific stereochemistry
    # [OH] connected to carbon at position 3 with beta stereochemistry
    # The [H] indicates explicit hydrogen, helps define stereochemistry
    # The @@ indicates beta configuration (below the plane)
    beta_3_oh_pattern = Chem.MolFromSmarts('[H][C@@]1[C@@H](O)CC[C@]2')
    
    # Alternative pattern that might match other valid 3beta-OH steroids
    alt_beta_3_oh_pattern = Chem.MolFromSmarts('[C@@H](O)CC[C@@]1')
    
    if not (mol.HasSubstructMatch(beta_3_oh_pattern) or mol.HasSubstructMatch(alt_beta_3_oh_pattern)):
        return False, "No 3beta-hydroxy group found"

    # Count rings to ensure we have a steroid-like structure
    ri = mol.GetRingInfo()
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check carbon count (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    return True, "Contains steroid core with 3beta-hydroxy group"