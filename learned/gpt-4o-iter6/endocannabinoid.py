"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for any long polyunsaturated chains with flexibility on length
    unsaturated_chain_pattern = Chem.MolFromSmarts("C=CCCC=CCCC=C")  # Partial structure
    long_chain_match = mol.HasSubstructMatch(unsaturated_chain_pattern)
    
    if not long_chain_match:
        return False, "No polyunsaturated long carbon chain typical of endocannabinoids found."
    
    # Check for key functional groups more comprehensively
    ethanolamine_pattern = Chem.MolFromSmarts("[NX3][CH2][CH2][OX2H]")
    amide_pattern = Chem.MolFromSmarts("N[C]=O")
    glycerol_pattern = Chem.MolFromSmarts("[OX2H][CX4]([OX2H])[CX4][OX2H]")
    has_key_group = mol.HasSubstructMatch(ethanolamine_pattern) or \
                    mol.HasSubstructMatch(amide_pattern) or \
                    mol.HasSubstructMatch(glycerol_pattern)

    if not has_key_group:
        return False, "No identified endocannabinoid functional group found (ethanolamine, amide, glycerol)."

    # Verify typical linkages more specifically
    ether_pattern = Chem.MolFromSmarts("COC")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    amide_linkage_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    
    has_linkages = mol.HasSubstructMatch(ether_pattern) or \
                   mol.HasSubstructMatch(ester_pattern) or \
                   mol.HasSubstructMatch(amide_linkage_pattern)

    if not has_linkages:
        return False, "No ether, ester, or amide linkages characteristic of endocannabinoids found."
    
    return True, "Matches characteristics of known endocannabinoids."