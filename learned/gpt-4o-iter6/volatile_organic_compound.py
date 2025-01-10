"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is assumed to have a boiling point less than or equal to 250 degreeC, often characterized by lower molecular weight.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is considered a volatile organic compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # For approximation, assume a compound with mol_wt less than 300 g/mol could be a VOC
    if mol_wt <= 300:
        return True, f"Molecular weight {mol_wt} suggests potential volatility"
    else:
        return False, f"Molecular weight {mol_wt} is too high for typical VOC classification"

# Test the function with a SMILES string
# Example: 1-dodecene (SMILES: CCCCCCCCCCC=C), expected to be a VOC
example_smiles = "CCCCCCCCCCC=C"
is_volatile, reason = is_volatile_organic_compound(example_smiles)
print(is_volatile, reason)