"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:38181 catechol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol component based on its SMILES string.
    A catechol component is defined as an o-diphenol group, which may be part of a larger ring system or have additional substituents.

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
    
    # Define a broader catechol pattern
    catechol_pattern = Chem.MolFromSmarts("[c;r3,r4,r5,r6]1:c(:c(:c(:c(:c:1)-[OH]):c-[OH])-*)-*")
    
    # Recursive function to find all potential catechol components
    def find_catechols(mol, catechol_matches):
        # Check if the molecule contains the catechol pattern
        matches = mol.GetSubstructMatches(catechol_pattern)
        catechol_matches.extend(matches)
        
        # Remove the identified catechol components from the molecule
        if matches:
            for match in matches:
                mol = Chem.DeleteSubstructs(mol, match)
            
            # Recursively search the remaining molecule
            find_catechols(mol, catechol_matches)
    
    catechol_matches = []
    find_catechols(mol, catechol_matches)
    
    if catechol_matches:
        # Additional checks or heuristics can be added here
        # For example, check if the catechol components are part of a larger ring system
        # or if they are connected to specific functional groups
        
        # Count the number of catechol components
        n_catechols = len(catechol_matches)
        
        # Check if the molecule has other non-catechol aromatic rings
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        non_catechol_rings = aromatic_rings - n_catechols
        
        # Heuristic: Molecules with multiple catechol components and no other aromatic rings are more likely to be true catechols
        if n_catechols > 1 and non_catechol_rings == 0:
            return True, f"Contains {n_catechols} catechol components"
        
        return True, "Contains a catechol component"
    
    return False, "Does not contain a catechol component"