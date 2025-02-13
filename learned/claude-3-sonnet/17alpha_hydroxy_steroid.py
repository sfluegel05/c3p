"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core pattern that allows for different bond types
    # and variations in the D ring
    steroid_core = Chem.MolFromSmarts(
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6](@[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~1)([#6,H])[#6]"
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 17-alpha-OH
    # Matches C17 with alpha OH and typical steroid connectivity
    oh_17_alpha = Chem.MolFromSmarts(
        '[C;R1]12[C;R1][C;R1][C;R1]3[C;R1]1[#6][#6][#6]1[#6][#6][#6][#6][C;R1]1[C@@]2([C;R1]3)[O;H1]'
    )
    
    # Alternative pattern for cases with explicit H
    oh_17_alpha_alt = Chem.MolFromSmarts(
        '[C;R1]12[C;R1][C;R1][C;R1]3[C;R1]1[#6][#6][#6]1[#6][#6][#6][#6][C;R1]1[C@]2([C;R1]3)[O;H1]'
    )

    if not (mol.HasSubstructMatch(oh_17_alpha) or mol.HasSubstructMatch(oh_17_alpha_alt)):
        return False, "No 17-alpha hydroxyl group found"

    # Basic validation checks
    # Count carbons (steroids typically have 19+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 19:
        return False, "Too few carbons for steroid structure"

    # Count rings (steroids have 4 fused rings)
    ri = mol.GetRingInfo()
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check for reasonable molecular weight (typical steroids are 250-500 Da)
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1000:  # Upper limit increased to accommodate glycosides
        return False, "Molecular weight outside typical range for steroids"

    # Count oxygens (should have at least the 17-OH oxygen)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found"

    # Verify sp3 carbons typical in steroid core
    sp3_carbons = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() == 6 and 
                     atom.GetHybridization() == Chem.HybridizationType.SP3)
    if sp3_carbons < 8:
        return False, "Insufficient sp3 carbons for steroid structure"

    return True, "Contains steroid core with 17-alpha hydroxyl group"