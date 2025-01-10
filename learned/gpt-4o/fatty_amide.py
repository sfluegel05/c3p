"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the amide group pattern: should have connectivity typical of fatty amides
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Check for contiguous long carbon chain typical of fatty acids
    long_chain_len = 12 # Consider chain length typical of a fatty acid
    carbon_chain_pattern = Chem.MolFromSmarts(f"[C]{{{long_chain_len},}}") # Dummy pattern
  
    # Perform substructure search for long carbon chains, allow flexible matching
    matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not matches:
        return False, "No sufficiently long carbon chain detected"

    # Now ensure that at least one chain is part of an amide motif
    for match in matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 7: # Nitrogen atom neighboring in amide bond
                    if mol.HasSubstructMatch(amide_pattern):
                        return True, "Contains an amide group and a long carbon chain characteristic of fatty amides"

    return False, "Amide group not correctly linked to a long carbon chain"