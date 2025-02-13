"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids typically feature a long hydrocarbon chain with multiple double bonds,
    often linked by amide or ester groups to small molecules like ethanolamine or glycerol.

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

    # Check for long hydrocarbon chains with at least 3 double bonds
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX2,CX1]~[CX4,CX3,CX2]~[CX2]!@[CX2]~[CX2]~[CX4,CX3,CX2,CX1]")
    chain_matches = mol.GetSubstructMatches(hydrocarbon_chain_pattern)
    if len(chain_matches) < 1:
        return False, f"Insufficient long hydrocarbon chains with multiple double bonds, found {len(chain_matches)}"

    # Check for amide linkage (common in ethanolamide such as in anandamide)
    amide_pattern = Chem.MolFromSmarts("NC=O")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Check for ester linkage (common in glycerol esters)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not amide_matches and not ester_matches:
        return False, "No linking groups found (amide or ester expected)"

    # Check for ethanolamine or glycerol fragments
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("OCCO")
    if not mol.HasSubstructMatch(ethanolamine_pattern) and not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No ethanolamine or glycerol fragments found"

    return True, "Possesses typical features of endocannabinoids including polyunsaturated chain and linkage to small molecules"

# Example testing of the function
test_smiles = "CCCCCCCCCCCCCCC(=O)NCCO"  # Example SMILES for palmitoyl ethanolamide
result, reason = is_endocannabinoid(test_smiles)
print(f"Classification: {result}, Reason: {reason}")