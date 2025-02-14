"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:28799 porphyrins
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_porphyrin(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins are defined as natural pigments containing a fundamental skeleton
    of four pyrrole nuclei united through the alpha-positions by four methine
    groups to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for porphyrin macrocycle pattern
    porphyrin_pattern = Chem.MolFromSmarts("[nH]1cccc1C=Cc1c[nH]c2c1cccc2C=Cc1c[nH]c2c1cccc2C=Cc1c[nH]c2c1cccc2C=C1")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin macrocycle found"
    
    # Check planarity
    planar = rdMolDescriptors.CalcPlanarityFast(mol)
    if planar < 0.9:
        return False, "Structure is not planar"
    
    # Check for optional metal coordination
    metal_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 20]
    if metal_atoms:
        if len(metal_atoms) > 1:
            return False, "Multiple metal atoms detected"
        metal_symbol = Chem.GetPeriodicTable().GetElementSymbol(metal_atoms[0])
        return True, f"Porphyrin macrocycle with coordinated {metal_symbol} atom"
    
    return True, "Porphyrin macrocycle"