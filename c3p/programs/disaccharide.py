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

    # Check for the presence of exactly two ring structures (suggests two sugar units)
    if rdMolDescriptors.CalcNumRings(mol) != 2:
        return False, "Does not contain exactly two ring structures"
        
    # A pattern suggesting a glycosidic bond: an ether linkage where one oxygen connects two sugar rings
    glycosidic_pattern = Chem.MolFromSmarts("[C;R1]1[O;X2]!@[C;R1]2")  # Smart pattern for acetal/ketal linkages between two rings
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond pattern found"

    # Count atoms to verify absence of non-sugar elements; typical disaccharides have carbon and oxygen with possible functional groups
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Assuming typical disaccharide mass range based on carbon and oxygen count
    if carbon_count < 10 or carbon_count > 24:
        return False, f"Unexpected number of carbons: {carbon_count} for disaccharide"
    if oxygen_count < 6 or oxygen_count > 12:
        return False, f"Unexpected number of oxygens: {oxygen_count} for disaccharide"

    # Check for sugars-like hydroxyl groups
    hydroxyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]"))
    if len(hydroxyl_groups) < 5:
        return False, f"Too few hydroxyl groups for a disaccharide structure, found: {len(hydroxyl_groups)}"
    
    return True, "Contains features characteristic of disaccharides: two sugar units linked by a glycosidic bond"