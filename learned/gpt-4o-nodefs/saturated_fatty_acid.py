"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a long aliphatic chain with a carboxylic acid group at one end
    and no unsaturation or rings.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group at terminal position
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Check for presence of any double bonds (unsaturation)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        return False, "Contains double bonds indicating unsaturation"
    
    # Check for triple bonds 
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C#C")):
        return False, "Contains triple bonds indicating unsaturation"
    
    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Rings detected, should be aliphatic"
    
    # Verify chain length meets minimum criteria for saturated fatty acids
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 8:
        return False, "Chain length too short for a saturated fatty acid"

    # Ensure absence of other functional groups or heteroatoms besides expected
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ["C", "O", "H"]:
            return False, "Unexpected heteroatoms found"

    return True, "Molecule is a saturated fatty acid"