"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogue
Definition: A molecule that can substitute for a normal nucleobase in nucleic acids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue must resemble a natural nucleobase (adenine, guanine, cytosine, thymine, uracil)
    but with modifications such as substitutions, additions, or alterations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for natural nucleobases
    adenine_pattern = Chem.MolFromSmarts("Nc1ncnc2ncnc12")
    guanine_pattern = Chem.MolFromSmarts("Nc1nc2c(=O)[nH]c(=O)nc2n1")
    cytosine_pattern = Chem.MolFromSmarts("Nc1nc(=O)[nH]cc1")
    thymine_pattern = Chem.MolFromSmarts("Cc1c(=O)[nH]c(=O)nc1")
    uracil_pattern = Chem.MolFromSmarts("O=c1cc[nH]c(=O)[nH]1")

    # Check if the molecule matches any natural nucleobase pattern
    if (mol.HasSubstructMatch(adenine_pattern) or
        mol.HasSubstructMatch(guanine_pattern) or
        mol.HasSubstructMatch(cytosine_pattern) or
        mol.HasSubstructMatch(thymine_pattern) or
        mol.HasSubstructMatch(uracil_pattern)):
        # If it matches, it is a nucleobase analogue
        return True, "Resembles a natural nucleobase with modifications"
    else:
        # If it doesn't match, check for common modifications
        # Look for nitrogen-containing heterocycles (common in nucleobases)
        nitrogen_heterocycle_pattern = Chem.MolFromSmarts("[n]")
        if not mol.HasSubstructMatch(nitrogen_heterocycle_pattern):
            return False, "No nitrogen-containing heterocycle found"

        # Check for carbonyl groups (common in nucleobases)
        carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
        if not mol.HasSubstructMatch(carbonyl_pattern):
            return False, "No carbonyl group found"

        # Check for amine groups (common in nucleobases)
        amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")
        if not mol.HasSubstructMatch(amine_pattern):
            return False, "No amine group found"

        # If it has these features but doesn't match natural nucleobases, it might still be an analogue
        return True, "Resembles a nucleobase analogue with modifications"

    return False, "Does not resemble a nucleobase analogue"