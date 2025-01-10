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
    # - Chromene ring system with ketone at position 4
    # - Hydroxy group at position 3
    # - Phenyl substituent at position 2
    flavonol_pattern = Chem.MolFromSmarts('[$(O=C1c2ccccc2OC(c3ccccc3)=C1O)]')
    
    # Alternative more general pattern to catch variations
    flavonol_pattern2 = Chem.MolFromSmarts('[$(Oc1c(=O)c2ccccc2oc1-c1ccccc1)]')
    
    if not (mol.HasSubstructMatch(flavonol_pattern) or mol.HasSubstructMatch(flavonol_pattern2)):
        return False, "Missing flavonol core structure"

    # Verify presence of ketone at position 4 and OH at position 3
    # This pattern specifically looks for the O=C-C(O) arrangement
    essential_groups = Chem.MolFromSmarts('[$(O=C1c2ccccc2O[C@H](O)C1)]')
    if not mol.HasSubstructMatch(essential_groups):
        return False, "Missing required ketone at position 4 or hydroxy at position 3"

    # Count basic ring systems (should have at least 2 aromatic rings)
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 2:
        return False, "Missing required ring systems"

    # Count oxygen atoms (should have at least 3: 3-OH, 4=O, and ring O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms for flavonol structure"

    # Additional check for correct connectivity
    connectivity_pattern = Chem.MolFromSmarts('[$(c1cc(O)c2c(c1)oc(-c1ccccc1)c(O)c2=O)]')
    if not mol.HasSubstructMatch(connectivity_pattern):
        # Try alternative pattern for substituted variants
        alt_pattern = Chem.MolFromSmarts('[$(c1c(O)c2oc(-c3ccccc3)c(O)c(=O)c2cc1O)]')
        if not mol.HasSubstructMatch(alt_pattern):
            return False, "Incorrect connectivity pattern for flavonol"

    return True, "Contains flavonol core structure with 3-hydroxy group and required features"