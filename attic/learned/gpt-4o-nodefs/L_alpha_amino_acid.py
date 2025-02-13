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

    # L-alpha amino acids generally have a chiral center (C) attached to NH2 and COOH
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[N][C@@H]([C])[C](=O)O")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "Does not have the L-alpha-amino acid framework"

    # Identify all chiral centers and ensure they match L-amino acid configuration
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    has_L_configuration = False
    
    for center in chiral_centers:
        atom_index, chirality = center
        atom = mol.GetAtomWithIdx(atom_index)
        
        # Check for specific C@@ signifying S-configuration equivalent to natural L
        if chirality == 'S' and atom.GetSymbol() == 'C':
            has_L_configuration = True
            break
    
    if not has_L_configuration:
        return False, "No chiral center with L-configuration (S) found"

    return True, "Contains L-alpha-amino acid structure with correct chirality"

# Example test
example_smiles = "N[C@@H](Cc1c(O)[nH]c2ccccc12)C(O)=O"
result, reason = is_L_alpha_amino_acid(example_smiles)
print(result, reason)