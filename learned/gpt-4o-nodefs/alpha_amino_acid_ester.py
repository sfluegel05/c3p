"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester has an alpha-carbon bonded to an amino group and an esterified carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for alpha-amino acid ester
    # - A carbon atom with a primary amine NH2 adjacent and an ester linkage
    amino_acid_ester_pattern = Chem.MolFromSmarts('[NX3][CX4H]([#6])[CX3](=O)[O][#6]')

    # Check if the molecule has the alpha-amino acid ester substructure
    if not mol.HasSubstructMatch(amino_acid_ester_pattern):
        return False, "No alpha-amino acid ester pattern found"

    # Check for the ester group specifically
    ester_group_pattern = Chem.MolFromSmarts('C(=O)O[#6]')
    ester_matches = mol.GetSubstructMatches(ester_group_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Check for the amino group being primary 
    primary_amino_pattern = Chem.MolFromSmarts('[NX3H2]')
    if not mol.HasSubstructMatch(primary_amino_pattern):
        return False, "No primary amino group found"

    # Successful match means the molecule fulfills the basic structure of an alpha-amino acid ester
    return True, "Molecule contains an alpha-amino acid ester structure"

# Example usage:
example_smiles = "COC(=O)CN"
result, reason = is_alpha_amino_acid_ester(example_smiles)
print(f"Classification result: {result}, Reason: {reason}")