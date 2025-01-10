"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid typically has an amino group, a carboxyl group, and a side chain
    attached to a central α-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for amino group and carboxyl group
    amino_pattern = Chem.MolFromSmarts("[NH2,NH1,NH0+]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    
    # Search for amino and carboxyl groups in the molecule
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not amino_matches:
        return False, "No amino group found"
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Check if amino and carboxyl groups are attached to the same carbon atom
    # Assume the patterns matched a single atom for amino and one or two for carboxyl (O only in case of salts)
    possible_a_carbon = {a[0] for a in amino_matches}.intersection({c[0] for c in carboxyl_matches})
    
    if not possible_a_carbon:
        return False, "Amino and carboxyl groups are not connected to the same α-carbon"
    
    # If the structure passes all tests, consider it an amino acid
    return True, "Identified as an amino acid with appropriate functional groups"