"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems are characterized by the presence of a bicyclic core structure
    consisting of a 4-membered beta-lactam ring fused to a 5-membered ring 
    containing sulfur, with typical stereochemistry and substitutions.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) True for carbapenem with reason, False with reason otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern for beta-lactam ring: 4-membered ring containing nitrogen and carbonyl group
    beta_lactam_pattern = Chem.MolFromSmarts("N1C(=O)C(C)C1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Improved pattern for 5-membered ring containing sulfur, fused to beta-lactam
    fused_ring_pattern = Chem.MolFromSmarts("C1SCC2N=C1C=O")
    if not mol.HasSubstructMatch(fused_ring_pattern):
        return False, "No sulfur-containing ring fused to beta-lactam found"

    return True, "Contains the characteristic structure of a carbapenem antibiotic"