"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is a fatty acid ester in which the carboxylic acid component is lauric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the lauroyl ester pattern (C12 chain with carbonyl and oxygen)
    # Including SMILES checks to see if a dodecanoic ester (-C(=O)O-) is present
    lauroyl_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")  # C12 chain
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O")

    # Check for the lauroyl chain and the ester linkage pattern
    if mol.HasSubstructMatch(lauroyl_chain_pattern) and mol.HasSubstructMatch(ester_linkage_pattern):
        # Double-check: Ensure ester linkage connects the C12 chain (context check)
        chains = mol.GetSubstructMatches(lauroyl_chain_pattern)
        esters = mol.GetSubstructMatches(ester_linkage_pattern)
        for chain in chains:
            for ester in esters:
                # Check if the ester oxygen is connected to the C12 chain.
                if mol.GetAtomWithIdx(ester[2]).GetIdx() == chain[-1]:  # Assert ester oxygen connectivity to C12 terminus
                    return True, "Contains characteristic lauroyl ester group indicative of dodecanoate ester"
                    
    return False, "Does not contain characteristic lauroyl ester group"