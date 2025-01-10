"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    
    A 2-oxo monocarboxylic acid should typically have a carbon (second carbon) doubly bonded to an oxygen 
    (ketone group) and connected to a carboxylic acid group (COOH).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined pattern for the 2-oxo group (RCOC(=O) as the principal component of the chain)
    oxo_pattern = Chem.MolFromSmarts("[#6][#6](=O)[#6]")  
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing 2-oxo group (ketone group at correct position)"
        
    # Refined pattern for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"

    # Additional check: ensure structural integrity of 2-oxo monocarboxylic framework
    for match in mol.GetSubstructMatches(oxo_pattern):
        oxo_carbon_idx = match[1]  # selects second carbon in the pattern
        for nbr in mol.GetAtomWithIdx(oxo_carbon_idx).GetNeighbors():
            if nbr.GetSymbol() == "C":  # ensure ketone carbon links to another carbon
                if any(neigh.GetSymbol() == "O" and neigh.GetNeighbors()[0].GetSymbol() == "C" for neigh in nbr.GetNeighbors()):
                    return True, "Contains 2-oxo group and carboxylic acid group"
    
    return False, "Could not match typical 2-oxo monocarboxylic acid framework"