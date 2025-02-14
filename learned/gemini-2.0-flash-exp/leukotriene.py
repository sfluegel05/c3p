"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is a C20 polyunsaturated fatty acid derivative with four double bonds,
    three of which are conjugated, derived from arachidonic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a leukotriene, False otherwise
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a 20-carbon chain (backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Less than 20 carbons; not a leukotriene"
    
    
    # 2. Check for four double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    if len(double_bond_matches) < 4:
         return False, f"Less than 4 double bonds; found {len(double_bond_matches)}"

    # 3. Check for a conjugated triene. 
    # This will catch the most common arrangement, but some leukotrienes can have this broken by a single bond.
    conjugated_triene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C") 
    triene_matches = mol.GetSubstructMatches(conjugated_triene_pattern)
    if not triene_matches:
         return False, "No conjugated triene (C=C-C=C-C=C) found"
    
    #4 Check for at least one hydroxyl and one carboxylic acid
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")

    if not mol.HasSubstructMatch(hydroxyl_pattern):
         return False, "No hydroxyl group found"
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Meets the criteria for a leukotriene: 20 carbons, four double bonds with 3 conjugated, a hydroxyl and a carboxylic group."