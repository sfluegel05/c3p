"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids feature long polyunsaturated hydrocarbon chains and are linked to small molecules like ethanolamine (amide), glycerol (ester), or similar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an endocannabinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for polyunsaturated long hydrocarbon chains (usually C20 with multiple cis-double bonds)
    long_chain_pattern = Chem.MolFromSmarts("[C]1(~[CH]=[CH]~[CH]=[CH]~[CH]=[CH]~[CH]=[CH]~CCCCC)~1")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No typical long, polyunsaturated chain found"

    # Check for amide linkage (common in ethanolamine such as in anandamide)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)C")
    if mol.HasSubstructMatch(amide_pattern):
        return True, "Possesses amide linkage typical of endocannabinoids"

    # Check for ester linkage (common in glycerol esters)
    ester_pattern = Chem.MolFromSmarts("[O][CX3](=O)C")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Possesses ester linkage typical of endocannabinoids"

    # Check for potential small alcohol linkages (ethanolamine, glycerol)
    small_molecule_linkage_pattern = Chem.MolFromSmarts("[O][CH2][CH2]O")
    if mol.HasSubstructMatch(small_molecule_linkage_pattern):
        return True, "Contains small molecule linkages typically found in endocannabinoids"

    return False, "Does not possess key structural features of endocannabinoids"

# Example testing of the function
test_smiles = "CCCCCCCCCCCCCCC(=O)NCCO"  # Example SMILES for palmitoyl ethanolamide
result, reason = is_endocannabinoid(test_smiles)
print(f"Classification: {result}, Reason: {reason}")