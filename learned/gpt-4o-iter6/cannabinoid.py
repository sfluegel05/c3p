"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoids
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    A cannabinoid typically contains long hydrocarbon chains, oxygen in a heterocyclic ring or 
    as part of functional groups, and often ester/amide/ether formations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for long hydrocarbon alkene/alkyne chains
    chain_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long hydrocarbon chain found"
    
    # Check for presence of oxygen function groups (hydroxyls, carbonyl, ether, ester linkage)
    oxygen_pattern = Chem.MolFromSmarts("[OX2,OX1]")
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No oxygen-containing functional group found"
    
    # Check for presence of heteroatoms, indicating non-simple hydrocarbon
    heteroatoms = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [6, 1]:  # Exclude carbon and hydrogen
            heteroatoms.add(atom.GetAtomicNum())
    
    if len(heteroatoms) == 0:
        return False, "Lacks heteroatoms, not a cannabinoid"

    # Optional: Handle specific functionalities, e.g., amide linkage
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_pattern):
        return True, "Classified as cannabinoid by presence of amide linkage"
    
    return True, "Classified as cannabinoid by general structural features"