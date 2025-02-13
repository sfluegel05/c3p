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
    # - Chromene-4-one system (benzopyran-4-one)
    # - Hydroxy group at position 3
    # - Phenyl substituent at position 2
    # The pattern allows for various substitutions on both rings
    flavonol_core = Chem.MolFromSmarts('[$(O=C1c2c(O)c(O)cc(O)c2OC(c3ccccc3)=C1O),'
                                      '$(O=C1c2c(O)cc(O)cc2OC(c3ccccc3)=C1O),'
                                      '$(O=C1c2cc(O)c(O)cc2OC(c3ccccc3)=C1O)]')
    
    # More general pattern that captures the essential chromene-4-one with 3-OH
    basic_pattern = Chem.MolFromSmarts('O=C1C=C(O)Oc2ccccc12')
    
    if not mol.HasSubstructMatch(flavonol_core) and not mol.HasSubstructMatch(basic_pattern):
        return False, "Missing basic flavonol structure"

    # Check for the essential 3-hydroxy and 4-keto arrangement
    # This pattern is more flexible than the previous version
    keto_hydroxy = Chem.MolFromSmarts('[$(O=C1c2c(-[*])cccc2OC(c3ccccc3)=C1O)]')
    if not mol.HasSubstructMatch(keto_hydroxy):
        return False, "Missing required 3-hydroxy and 4-keto groups"

    # Verify presence of two ring systems (AC rings of flavonoid skeleton)
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 2:
        return False, "Missing required ring systems"

    # Count oxygen atoms (should have at least 3: 3-OH, 4=O, and ring O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms for flavonol structure"

    # Additional check for aromatic character of rings
    aromatic_rings = Chem.MolFromSmarts('[$(c1ccccc1):1].[$(c1ccccc1):2]')
    if not mol.HasSubstructMatch(aromatic_rings):
        return False, "Missing required aromatic rings"

    return True, "Contains flavonol core structure with required 3-hydroxy group and ketone at position 4"