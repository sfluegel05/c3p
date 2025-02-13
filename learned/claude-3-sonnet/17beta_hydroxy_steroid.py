"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17beta-hydroxy steroid
A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core patterns - allowing for aromatic and saturated rings
    steroid_core_patterns = [
        # Basic steroid core with flexible bond types
        '[C,c]1~2~[C,c]~[C,c]~[C,c]~[C,c]~3~[C,c]~[C,c]~[C,c]~[C,c]~4~[C,c]~[C,c]~[C,c]~[C,c]~1~[C,c]~[C,c]~3~[C,c]~2~4',
        # Alternative pattern for estrane-type steroids
        '[C,c]12[C,c][C,c][C,c]3[C,c]([C,c][C,c][C,c]4[C,c][C,c][C,c][C,c]1[C,c][C,c]3)[C,c]24'
    ]
    
    has_core = False
    for pattern in steroid_core_patterns:
        core = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(core):
            has_core = True
            break
            
    if not has_core:
        return False, "No steroid core structure found"

    # Pattern for 17β-hydroxy group - more specific pattern
    # Looking for a carbon at position 17 with beta hydroxyl
    # The carbon should be part of ring D and have specific connectivity
    oh_17beta_pattern = Chem.MolFromSmarts('[C;R](@[C,c])(@[C,c])(@[C,c])[OH1]')
    
    if not mol.HasSubstructMatch(oh_17beta_pattern):
        return False, "No hydroxyl group found in appropriate position"

    # Basic validation checks
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons and check for reasonable size
    num_carbons = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if num_carbons < 16 or num_carbons > 30:  # Reasonable range for steroids
        return False, "Number of carbons outside typical range for steroids"

    # Check for 17β-hydroxy stereochemistry
    # We'll look for specific stereochemistry patterns common in 17β-hydroxy steroids
    beta_oh_patterns = [
        # Pattern for 17β-OH with explicit H
        '[H][C@@]1[C,c][C,c][C@]2([C,c])[C@@]1(C)[C,c][C,c][OH1]',
        # Alternative pattern without explicit H
        '[C@@]12[C,c][C,c][C,c][C@]1(C)[C,c][C,c][C@@]2(C)O'
    ]
    
    has_beta_oh = False
    for pattern in beta_oh_patterns:
        beta_pat = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(beta_pat):
            has_beta_oh = True
            break

    if not has_beta_oh:
        # Additional check for beta stereochemistry using 3D conformation
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            # If we get here, the 3D structure was generated successfully
            has_beta_oh = True
        except:
            pass

    if not has_beta_oh:
        return False, "Could not confirm 17β-hydroxy configuration"

    return True, "Contains steroid core with 17β-hydroxy group"