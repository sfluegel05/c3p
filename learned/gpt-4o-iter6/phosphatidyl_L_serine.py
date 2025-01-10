"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Classifies a molecule as a phosphatidyl-L-serine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """

    # Convert the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    
    # Return false with a reason if the molecule object cannot be created.
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the phosphatidyl fragment with bound serine.
    phosphatidyl_pattern = Chem.MolFromSmarts(
        "COP(=O)(O)O[C@H](COC(=O)C)C(OC(=O)C)O"
    )
    
    # Check for the presence of a phosphatidyl group in the molecule.
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl linkage with serine found"
    
    # Define the SMARTS pattern for the esterified serine (serine's hydroxy esterified).
    serine_pattern = Chem.MolFromSmarts("[N+][C@@H](CC(=O)O)C(=O)O")
    
    # Check for the presence of the serine esterified group.
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "Missing or unestablished serine group"
    
    # Define the SMARTS pattern for ester linkages indicative of fatty acid chains.
    ester_linkage_pattern = Chem.MolFromSmarts("OC(=O)CC")
    
    # Check if the molecule has at least two ester linkages.
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester linkages; found {len(ester_matches)}"

    return True, "Molecule meets the criteria for phosphatidyl-L-serine"