"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide with a phytosphingosine backbone 
    having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Phytosphingosine backbone improved pattern (stereocenters and hydroxyl groups)
    phytosphingosine_pattern = Chem.MolFromSmarts("C[C@H](O)[C@H](O)[C@@H](CCC(O)O)N")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"

    # Fatty acyl amide linkage pattern - long carbon chain with amide bond
    # at least 9 additional carbons in the chain (arbitrary threshold)
    acyl_amide_pattern = Chem.MolFromSmarts("C(=O)[NH1][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(acyl_amide_pattern):
        return False, "No long-chain fatty acyl amide group found attached to nitrogen"
    
    # Check the chain length (optional more stringent checks)
    chain_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)  # carbon count
    if chain_count < 20:  # arbitrary carbon threshold for N-acylphytosphingosines
        return False, "Insufficient carbon count for fatty acid chain"

    return True, "Molecule contains a phytosphingosine backbone with a fatty acyl group attached to nitrogen"