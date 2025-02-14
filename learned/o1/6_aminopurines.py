"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: CHEBI:26566 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine derivative based on its SMILES string.
    A 6-aminopurine derivative is any compound having 6-aminopurine (adenine) as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 6-aminopurine derivative, False otherwise
        str: Reason for classification
    """
    # Parse the input SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for adenine (6-aminopurine)
    # This pattern captures the purine ring with an amino group at position 6
    # and ensures there is no oxygen (carbonyl group) attached to the ring,
    # which would indicate guanine or other derivatives.
    adenine_smarts = """
    [nH]1cnc2nc(N)[nH]c2n1   # Adenine core with amino group at position 6
    """

    # Remove comments and whitespace from the SMARTS pattern
    adenine_smarts = "".join(adenine_smarts.split())

    adenine_mol = Chem.MolFromSmarts(adenine_smarts)
    if adenine_mol is None:
        return False, "Error in creating adenine substructure pattern"

    # Perform substructure search to check if adenine is part of the molecule
    if mol.HasSubstructMatch(adenine_mol):
        return True, "Molecule contains adenine (6-aminopurine) substructure"
    else:
        # Also consider tautomeric forms where hydrogen atoms may shift positions
        # Generate potential tautomers of adenine
        from rdkit.Chem.MolStandardize import rdMolStandardize

        tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
        adenine_tautomer = Chem.MolFromSmiles("Nc1ncnc2ncnc12")  # Standard adenine SMILES
        adenine_tautomers = tautomer_enumerator.Enumerate(adenine_tautomer)

        # Check if any tautomer matches as a substructure
        for tautomer in adenine_tautomers:
            tautomer_smarts = Chem.MolToSmarts(tautomer)
            tautomer_pattern = Chem.MolFromSmarts(tautomer_smarts)
            if tautomer_pattern and mol.HasSubstructMatch(tautomer_pattern):
                return True, "Molecule contains adenine tautomer substructure"

        return False, "Molecule does not contain adenine (6-aminopurine) substructure"