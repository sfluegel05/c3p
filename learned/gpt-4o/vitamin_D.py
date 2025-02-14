"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    
    Vitamins D are identified by their secosteroid backbone and specific side chains
    and hydroxyl groups typical for this class of compounds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broader secosteroid recognition pattern: open B-ring as in vitamin D
    # Pattern represents a sterane backbone with an open B-ring
    seco_pattern = Chem.MolFromSmarts("C1C[C@H](CCC2=CC=C3[C@@H]2CCCC3=C)CC1")
    if not mol.HasSubstructMatch(seco_pattern):
        return False, "No characteristic secosteroid structure found"
    
    # Check for the presence of multiple hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)C")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_count < 2:
        return False, f"Insufficient hydroxyl group matches, found {hydroxyl_count}"
    
    # Check for long alkyl chain that is halted at C=C or alcohols, typical in Vitamin D
    alkyl_chain_pattern = Chem.MolFromSmarts("CCCCCC(O)=C")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "Missing typical alkyl side chain termination"

    return True, "Molecule matches key structural features of vitamin D"