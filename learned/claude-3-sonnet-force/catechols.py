"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:38181 catechol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol component based on its SMILES string.
    A catechol component is defined as an o-diphenol group, without additional substituents on the benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a catechol component, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more specific catechol pattern
    catechol_pattern = Chem.MolFromSmarts("[c;r5]1:c(:c(:c(:c(:c:1)-[OH]):c-[OH])-*)-*")
    
    # Check if the molecule contains the catechol pattern
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    if catechol_matches:
        # Additional checks to ensure a valid catechol component
        for match in catechol_matches:
            # Check if the catechol component is part of a fused ring system
            if not mol.GetAtomWithIdx(match[0]).IsInRingOfSize(6):
                continue
            
            # Check if the catechol component has additional substituents on the benzene ring
            ring_atoms = mol.GetAtomRingData(match[0])
            if any(atom.GetTotalNumHs() < 1 for atom in ring_atoms):
                continue
            
            # Additional checks or heuristics can be added here
            
            # If all checks pass, classify as a catechol
            return True, "Contains a catechol component"
    
    return False, "Does not contain a catechol component"