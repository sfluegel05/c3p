"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the basic structure
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)[C](O)=O")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "Does not have the L-alpha-amino acid framework"

    # Check for chirality (L-configuration)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not any(center[1] == 'R' for center in chiral_centers):
        return False, "No chiral center with L-configuration (R) found"

    return True, "Contains L-alpha-amino acid structure with correct chirality"

# Example test
example_smiles = "N[C@@H](Cc1c(O)[nH]c2ccccc12)C(O)=O"
result, reason = is_L_alpha_amino_acid(example_smiles)
print(result, reason)