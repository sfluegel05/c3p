"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Glucosinolates have a thioglucoside linkage, a central C linked via S
    and N to a sulfonated oxime group, and carry a side-group.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a glucosinolate, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a thioglycoside structure: sugar ring linked through sulfur
    thioglycoside_pattern = Chem.MolFromSmarts("S[C@H]1O[C@@H]([C@H](O)[C@H](O)[C@@H]1O)CO")
    if not mol.HasSubstructMatch(thioglycoside_pattern):
        return False, "No thioglycoside linkage (linked through sulfur) found"

    # Check for sulfonated oxime group linked to central carbon: N=C-OS(=O)(=O)[O-]
    oxime_sulfonate_pattern = Chem.MolFromSmarts("N=C/[CX3](=[NX2]OS(=O)(=O)[O-])")
    if not mol.HasSubstructMatch(oxime_sulfonate_pattern):
        return False, "No appropriate sulfonated oxime linkage found"

    # Check for a sulfur linkage to central carbon: C-S-N and sidechain
    central_c_pattern = Chem.MolFromSmarts("[#6]-[S]-[#6](=NOS(=O)(=O)[O-])[CX4,CX3]")
    if not mol.HasSubstructMatch(central_c_pattern):
        return False, "Central carbon structure not recognized with correct linkages"
    
    return True, "Contains all structural elements of a glucosinolate"