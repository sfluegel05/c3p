"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: CHEBI:25394 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine contains an amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the monoamine pattern: aromatic ring -> 2-carbon chain -> amino group
    monoamine_pattern = Chem.MolFromSmarts("[c]CC[NH2,NH,N;!$(NC=O);!$(NS(=O)=O);!$(N[C,S]=O)]")
    
    # Check if the pattern matches
    if not mol.HasSubstructMatch(monoamine_pattern):
        return False, "No monoamine pattern found (aromatic ring -> 2C chain -> amino group)"

    # Verify the aromatic ring is present and not part of a larger functional group
    aromatic_pattern = Chem.MolFromSmarts("[c]")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"

    # Verify the amino group is present and not part of a larger functional group
    amino_pattern = Chem.MolFromSmarts("[NH2,NH,N;!$(NC=O);!$(NS(=O)=O);!$(N[C,S]=O)]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # Verify the two-carbon chain is present
    carbon_chain_pattern = Chem.MolFromSmarts("CC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No two-carbon chain found"

    # Additional check to ensure the amino group is not part of a larger functional group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            if atom.GetDegree() > 1:  # If nitrogen is connected to more than one atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:  # Carbon neighbor
                        if neighbor.GetDegree() > 2:  # If carbon is part of a larger functional group
                            return False, "Amino group is part of a larger functional group"

    return True, "Contains aromatic ring connected to amino group via two-carbon chain"