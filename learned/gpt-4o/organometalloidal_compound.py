"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determine if a molecule is an organometalloidal compound based on its SMILES string.
    Organometalloidal compounds have bonds between metalloids (e.g., As) and organic carbon atoms.

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
    
    # Arsenic directly bonded to carbon
    arsenic_carbon_pattern = Chem.MolFromSmarts("[As][C]")
    # Arsenic surrounded by carbon and oxygen, covering carboxylic/functional groups
    organoarsenic_pattern = Chem.MolFromSmarts("[C][As](=O)")
    # Arsenic within aromatic rings or complex structures
    aromatic_arsenic_pattern = Chem.MolFromSmarts("c[As]")

    # Checking for multiple arsenic atoms, potential polymeric features
    multiple_arsenic_pattern = Chem.MolFromSmarts("[As][As]")
    
    found_patterns = []
    
    # Check for each pattern and record reason
    if mol.HasSubstructMatch(arsenic_carbon_pattern):
        found_patterns.append("Contains arsenic-carbon bonds")
    if mol.HasSubstructMatch(organoarsenic_pattern):
        found_patterns.append("Contains functional organyl-arsenic component")
    if mol.HasSubstructMatch(aromatic_arsenic_pattern):
        found_patterns.append("Contains arsenic in an aromatic context")
    if mol.HasSubstructMatch(multiple_arsenic_pattern):
        found_patterns.append("Contains multiple arsenic atoms forming complex bonds")

    if found_patterns:
        return True, " and ".join(found_patterns) + " indicative of organometalloidal compounds"

    return False, "Does not meet organometalloidal structural criteria"