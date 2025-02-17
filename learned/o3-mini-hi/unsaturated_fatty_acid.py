"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: Unsaturated Fatty Acid
Definition:
    Any fatty acid containing at least one C=C or C#C bond.
A fatty acid is assumed to have a carboxylic acid (COOH) group.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    A molecule is considered an unsaturated fatty acid if it includes:
      - A carboxylic acid functionality (e.g., COOH group)
      - At least one unsaturation in the form of a carbon-carbon double (C=C)
        or triple bond (C#C).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an unsaturated fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for a carboxylic acid group.
    # We consider both protonated (COOH) and deprotonated (COO-) forms.
    acid_pattern_prot = Chem.MolFromSmarts("C(=O)[O;H1]")    # protonated acid: C(=O)OH
    acid_pattern_deprot = Chem.MolFromSmarts("C(=O)[O-]")      # deprotonated acid: C(=O)O-

    has_acid = mol.HasSubstructMatch(acid_pattern_prot) or mol.HasSubstructMatch(acid_pattern_deprot)
    if not has_acid:
        return False, "Does not contain a carboxylic acid group"

    # Define SMARTS patterns for unsaturation:
    double_bond_pattern = Chem.MolFromSmarts("C=C")  # carbon-carbon double bond
    triple_bond_pattern = Chem.MolFromSmarts("C#C")   # carbon-carbon triple bond

    has_double = mol.HasSubstructMatch(double_bond_pattern)
    has_triple = mol.HasSubstructMatch(triple_bond_pattern)
    if not (has_double or has_triple):
        return False, "Contains no carbon–carbon double or triple bonds"

    return True, "Contains a carboxyl group and at least one C=C or C#C bond (unsaturation)"
    
# Example usage (uncomment below lines to test):
# test_smiles = "OC(=O)CCC=C"  # pent-4-enoic acid
# result, reason = is_unsaturated_fatty_acid(test_smiles)
# print(result, reason)