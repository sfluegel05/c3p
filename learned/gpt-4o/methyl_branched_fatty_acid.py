"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.

    A methyl-branched fatty acid is defined as any branched-chain fatty acid
    containing methyl branches only.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # Accommodating carboxylate anion as well
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Get all carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Check each carbon atom for methyl branches and ensure no other types of branches exist
    for carbon in carbon_atoms:
        num_hydrogens = sum(1 for neighbor in carbon.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        num_carbon_neighbors = sum(1 for neighbor in carbon.GetNeighbors() if neighbor.GetAtomicNum() == 6)
        
        # Check for methyl groups (CH3)
        if num_carbon_neighbors > 1:  # Branching occurs here
            branch_atom = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != carbon.GetIdx()]
            # Check if branching is exclusively via methyl groups
            for b in branch_atom:
                if sum(1 for n in b.GetNeighbors() if n.GetAtomicNum() == 1) != 3:
                    return False, f"Non-methyl branch found at carbon index {b.GetIdx()}"
    
    return True, "Molecule contains methyl branches and a carboxylic acid group, no other types of branches"