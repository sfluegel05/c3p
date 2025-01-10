"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Check for carboxylic acid functional group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # Consider deprotonated state as well
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Check carbon count to ensure it's within a reasonable range for fatty acids
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 4 or len(carbon_atoms) > 50:
        return False, "Unlikely carbon count for a fatty acid"
   
    # Define SMARTS patterns for a variety of branching patterns
    branching_patterns = [
        Chem.MolFromSmarts("[C](C)(C)"),  # Methyl branching
        Chem.MolFromSmarts("[C](C)(C)C"),  # More complex branches
        Chem.MolFromSmarts("C(C)(C)C"),   # Isopropyl
        Chem.MolFromSmarts("[C]=C"),  # Unsaturated branches
    ]

    # Check branching pattern matches to confirm it's a branched structure
    for pattern in branching_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains terminal carboxylic acid group and appropriate branching pattern"
    
    return False, "No suitable alkyl branching pattern found"