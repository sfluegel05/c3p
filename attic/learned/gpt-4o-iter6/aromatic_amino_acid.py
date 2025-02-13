"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid has both an aromatic ring and an amino acid structure (amino and carboxylate group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for an aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("a")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic ring found"
    
    # Look for amino acid structure (amino group NH2 and carboxylate C(=O)O)
    amino_group_pattern = Chem.MolFromSmarts("[NX3][CX4]")
    carboxylate_group_pattern = Chem.MolFromSmarts("[CX3](=O)[O]")
    
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No amino group found"
        
    if not mol.HasSubstructMatch(carboxylate_group_pattern):
        return False, "No carboxylic acid group found"

    # Check connectivity between aromatic ring and amino acid moiety
    # Aromatic ring connected to a carbon which is bonded to amino and carboxylate groups
    aromatic_amino_acid_pattern = Chem.MolFromSmarts("a-[CX4](N)[CX3](=O)O")
    if not mol.HasSubstructMatch(aromatic_amino_acid_pattern):
        return False, "Aromatic ring not connected to amino acid moiety"

    return True, "Contains aromatic ring and amino acid moiety connected appropriately"

# Test the function with a sample SMILES string
print(is_aromatic_amino_acid("NC(Cc1ccc(O)c(O)c1)C(O)=O"))  # Example for L-dopa