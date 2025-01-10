"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid zwitterion based on its SMILES string.

    An alpha-amino acid zwitterion is characterized by having a protonated amino group (NH3+)
    and a deprotonated carboxylic acid group (C(=O)[O-]).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for protonated amino group [NH3+]
    amino_pattern = Chem.MolFromSmarts("[NH3+]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No protonated amino group [NH3+] found"
    
    # Look for deprotonated carboxylate group C(=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No deprotonated carboxylate group C(=O)[O-] found"
    
    return True, "Molecule is an alpha-amino acid zwitterion with [NH3+] and C(=O)[O-]"

# Example usage:
# smiles = "[NH3+][C@@H](CCC(=O)NCCc1ccc(O)cc1)C([O-])=O"
# print(is_alpha_amino_acid_zwitterion(smiles))