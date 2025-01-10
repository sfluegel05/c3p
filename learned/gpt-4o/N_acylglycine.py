"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is characterized by an acyl group bonded as an amide to the nitrogen of a glycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycine pattern (NCC(=O)O) - represents the core glycine structure
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine core structure found"
    
    # Look for an amide linkage pattern (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No acyl group attached as amide linkage"

    # Get the substructure matches for both patterns
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Check if the nitrogen in glycine is part of an amide bond
    for amide_match in amide_matches:
        for glycine_match in glycine_matches:
            # Check if the glycine nitrogen is involved in an amide linkage
            if amide_match[1] == glycine_match[0]:  # amide N is the glycine N
                return True, "Found N-acylglycine characteristic structure"

    return False, "No acyl group bonded to glycine nitrogen as expected for N-acylglycine"