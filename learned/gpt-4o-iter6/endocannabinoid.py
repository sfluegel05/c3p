"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

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
    
    # Adjust criteria for long carbon chain identifying both straight and kinked structures
    polyunsaturated_chain_pattern = Chem.MolFromSmarts("C=CCCC=CCCC=CCC")
    long_chain_match = mol.HasSubstructMatch(polyunsaturated_chain_pattern)

    if not long_chain_match:
        return False, "No polyunsaturated long carbon chain typical of endocannabinoids found."
    
    # Check for specific key functional groups like ethanolamine, distinct patterns
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    amide_pattern = Chem.MolFromSmarts("NC=O")
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    has_key_group = mol.HasSubstructMatch(ethanolamine_pattern) or \
                    mol.HasSubstructMatch(amide_pattern) or \
                    mol.HasSubstructMatch(glycerol_pattern)

    if not has_key_group:
        return False, "No identified endocannabinoid functional group (ethanolamine, amide, glycerol) found."

    # Check for typical linkages
    ether_pattern = Chem.MolFromSmarts("COC")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_linkage_pattern = Chem.MolFromSmarts("C(=O)N")
    has_linkages = mol.HasSubstructMatch(ether_pattern) or \
                   mol.HasSubstructMatch(ester_pattern) or \
                   mol.HasSubstructMatch(amide_linkage_pattern)

    if not has_linkages:
        return False, "No ether, ester, or amide linkages characteristic of endocannabinoids found."
    
    return True, "Matches characteristics of known endocannabinoids."