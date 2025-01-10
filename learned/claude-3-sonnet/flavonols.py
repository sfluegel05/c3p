"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:28413 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    Flavonols are hydroxyflavones with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavonol core pattern:
    # - Chromene-4-one system with 2-phenyl and 3-hydroxy
    # [#6] represents any carbon (aliphatic or aromatic)
    # The pattern allows for various substitutions on both rings
    flavonol_core = Chem.MolFromSmarts('[$(O=C1[#6]2[#6][#6][#6][#6][#6]2O[#6](c3[#6][#6][#6][#6][#6]3)=C1O)]')
    
    if not mol.HasSubstructMatch(flavonol_core):
        return False, "Missing basic flavonol structure"

    # More specific check for the 3-hydroxy and 4-keto arrangement
    # This ensures the oxygen is at position 3 and the ketone at position 4
    keto_hydroxy = Chem.MolFromSmarts('[$(O=C1c2[#6][#6][#6][#6][#6]2OC(c3[#6][#6][#6][#6][#6]3)=C1O)]')
    if not mol.HasSubstructMatch(keto_hydroxy):
        return False, "Missing required 3-hydroxy and 4-keto groups in correct positions"

    # Verify the presence of required ring systems
    # Need benzopyran system (fused rings) and a phenyl ring
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 3:  # Need at least 3 rings (A, C, and B rings)
        return False, "Missing required ring systems"

    # Count oxygen atoms (should have at least 3: 3-OH, 4=O, and ring O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms for flavonol structure"

    # Additional check for aromatic character
    aromatic_pattern = Chem.MolFromSmarts('c1ccccc1')  # At least one benzene ring
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "Missing required aromatic system"

    # Check that the structure has the correct number of double bonds
    # This helps confirm the chromene-4-one system
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    if len(mol.GetSubstructMatches(double_bond_pattern)) < 1:
        return False, "Missing required double bonds"

    return True, "Contains flavonol core structure with required 3-hydroxy group and ketone at position 4"