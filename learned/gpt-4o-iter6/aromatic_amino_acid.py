"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid has both an aromatic ring and an amino acid structure (amino and carboxyl group).

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
    aromatic_ring_pattern = Chem.MolFromSmarts("[a]")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic ring found"
    
    # Look for amino group, specifically targeting connection that fits in overall amino acid structure
    amino_group_pattern = Chem.MolFromSmarts("[NX3][CX4]")
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No suitable amino group found"
    
    # Look for carboxyl group (C(=O)O)
    carboxyl_group_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxyl_group_pattern):
        return False, "No carboxyl group found"
    
    # Check the connectivity: The aromatic ring directly connects to a carbon atom that also connects to an amine and a carboxylate structure within a short chain (three bonds)
    aromatic_amino_acid_pattern = Chem.MolFromSmarts("a-[CX3]-[NX3]-C(=O)[O]")
    if not mol.HasSubstructMatch(aromatic_amino_acid_pattern):
        return False, "Aromatic ring not properly connected to amino acid functionality"

    return True, "Contains aromatic ring and amino acid moiety connected appropriately"

# Test the function with an example SMILES string for L-dopa
print(is_aromatic_amino_acid("NC(Cc1ccc(O)c(O)c1)C(O)=O"))