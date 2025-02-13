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

    # Check for long polyunsaturated chains typically seen in fatty acids
    polyunsaturated_chain_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(polyunsaturated_chain_pattern):
        has_long_chain = True
    else:
        return False, "No long polyunsaturated hydrocarbon chain found"

    # Check for amide linkage (common in ethanolamine such as in anandamide)
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    has_amide = mol.HasSubstructMatch(amide_pattern)

    # Check for ester linkage (common in glycerol esters)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    has_ester = mol.HasSubstructMatch(ester_pattern)

    # Check for small molecule linkages like ethanolamine or glycerol
    ethanolamine_linkage = Chem.MolFromSmarts("NCCO")
    glycerol_linkage = Chem.MolFromSmarts("OCCO")
    has_small_linkage = mol.HasSubstructMatch(ethanolamine_linkage) or mol.HasSubstructMatch(glycerol_linkage)

    if has_amide or has_ester or has_small_linkage:
        return True, "Contains structural features typical of endocannabinoids"

    return False, "Does not possess key structural features of endocannabinoids"

# Example testing of the function
test_smiles = "CCCCCCCCCCCCCCC(=O)NCCO"  # Example SMILES for palmitoyl ethanolamide
result, reason = is_endocannabinoid(test_smiles)
print(f"Classification: {result}, Reason: {reason}")