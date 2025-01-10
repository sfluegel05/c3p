"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid contains an aromatic ring connected to an amino and carboxyl group.

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
    
    # Correct pattern to find the primary amine connected to a carbon 
    amino_group_pattern = Chem.MolFromSmarts("[NX3;H2]-[CX4][CX3]=O")
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No suitable amino group found"
    
    # Confirm the presence of a carboxyl group
    carboxyl_group_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_group_pattern):
        return False, "No carboxyl group found"

    # Check for proper connectivity (aromatic ring attached to the amino acid backbone)
    aromatic_amino_acid_pattern = Chem.MolFromSmarts("a-[CH2]-[NH2]-C(=O)O")
    if not mol.HasSubstructMatch(aromatic_amino_acid_pattern):
        return False, "Aromatic ring not properly connected to amino acid functionality"
   
    return True, "Contains aromatic ring and amino acid moiety connected appropriately"

# Test the function with an example SMILES string
print(is_aromatic_amino_acid("NC(Cc1ccc(O)c(O)c1)C(O)=O"))