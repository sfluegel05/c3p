"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is defined as two monosaccharides joined by a glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adjusted boundaries for the number of rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 1 or num_rings > 3:
        return False, f"Number of rings ({num_rings}) may signal complexity beyond disaccharides, need 1-3"

    # Check for the presence of glycosidic bonds broadly
    glycosidic_patterns = [
        Chem.MolFromSmarts("O[C;R][C;R]O"),  # classic O-linked
        Chem.MolFromSmarts("O[C;!R][!R]C"),  # alternative linkage
        Chem.MolFromSmarts("COC"),  # Simple ether backbone
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns):
        return False, "No clear glycosidic bond pattern found"

    # Count of carbon and oxygen atoms expanded
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if carbon_count < 10 or carbon_count > 36:
        return False, f"Unexpected number of carbons: {carbon_count}, typical for disaccharides range 10-36"
    if oxygen_count < 5 or oxygen_count > 18:
        return False, f"Unexpected number of oxygens: {oxygen_count}, typical for disaccharides range 5-18"

    # Checking for hydroxyl groups to confirm sugar structure
    num_hydroxyls = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]")))
    if num_hydroxyls < 4:
        return False, f"Too few hydroxyl groups ({num_hydroxyls}) for a typical sugar structure"

    return True, "Contains features characteristic of disaccharides: two sugar units linked by a glycosidic bond"