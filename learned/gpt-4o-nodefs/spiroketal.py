"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a given SMILES string represents a spiroketal.
    
    A spiroketal is characterized by:
    - A spiro junction
    - At least one cyclic ketal group
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a spiroketal, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a spiro junction: two rings sharing a common atom
    spiro_pattern = Chem.MolFromSmarts('[R]@[R]')
    if not mol.HasSubstructMatch(spiro_pattern):
        return False, "No spiro junction found"
    
    # Define the SMARTS pattern for a ketal group: C-O-C-O-C
    ketal_pattern = Chem.MolFromSmarts('[C]([O])[C]([O])[C]')
    if not mol.HasSubstructMatch(ketal_pattern):
        return False, "No ketal group found"
    
    # Verify the presence of cyclic (part of a ring) structure in ketal
    ring_info = mol.GetRingInfo()
    ring_bonds = ring_info.BondRings()
    ketal_matches = mol.GetSubstructMatches(ketal_pattern)
    
    if not any(set(m) <= set(bond) for bond in ring_bonds for m in ketal_matches):
        return False, "No cyclic ketal group found"
    
    return True, "Contains a spiro junction and a cyclic ketal group"