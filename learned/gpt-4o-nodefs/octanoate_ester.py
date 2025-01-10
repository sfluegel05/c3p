"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester contains an ester linkage where the alcohol moiety
    has formed an ester with octanoic acid, characterized by an 8-carbon aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an ester linkage pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check for basic ester bond
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond found"

    # Find all substructures that match the ester pattern
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    for match in ester_matches:
        # Get the carboxylate carbon atom index
        carboxylate_carbon_idx = match[0]
        
        # Check the aliphatic chain starting from the carboxylate carbon
        carbon_chain = Chem.rdmolops.GetShortestPath(mol, carboxylate_carbon_idx, match[1])
        
        # If there are exactly 8 carbon atoms in the aliphatic chain, it's an octanoate
        if len([atom for atom in carbon_chain if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6]) == 8:
            return True, "Molecule contains an octanoate ester group"

    return False, "No octanoate ester group with 8-carbon chain found"