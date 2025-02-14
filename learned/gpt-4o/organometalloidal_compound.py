"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determine if a molecule is an organometalloidal compound based on its SMILES string.
    Organometalloidal compounds have bonds between a metalloid (e.g., As) and organic carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for arsenic-carbon bonds
    arsenic_carbon_pattern = Chem.MolFromSmarts("[As]-[C]")
    if mol.HasSubstructMatch(arsenic_carbon_pattern):
        # Confirm the presence of arsenic and carbon bonds
        return True, "Contains arsenic-carbon bonds indicative of organometalloidal compounds"
    
    return False, "No arsenic-carbon bonds found"