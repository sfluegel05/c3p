"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid-zwitterion obtained by transfer of a proton from the carboxy to the amino group of any alpha-amino acid; major species at pH 7.3.
The program checks if the molecule contains at least one alpha carbon bearing both a protonated amine ([NH3+]) and a deprotonated carboxylate (C(=O)[O-]) group.
"""

from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    
    The method looks for the characteristic substructure in which an alpha carbon is attached
    both to a protonated amino group ([NH3+]) and to a deprotonated carboxylate group (C(=O)[O-]).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an alpha-amino-acid zwitterion, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern representing an alpha carbon attached to [NH3+] and C(=O)[O-]
    # This pattern will match a carbon (irrespective of chirality marks) that is bound to a protonated amine and a carboxylate group.
    # It covers the key features of an alpha-amino-acid zwitterion.
    aa_zwitterion_pattern = Chem.MolFromSmarts("[C]([NH3+])C(=O)[O-]")
    
    if aa_zwitterion_pattern is None:
        return False, "Failed to create SMARTS pattern for alpha-amino-acid zwitterion"
    
    # Search for the substructure match.
    if not mol.HasSubstructMatch(aa_zwitterion_pattern):
        return False, "The key substructure [C]([NH3+])C(=O)[O-] representing an alpha-amino-acid zwitterion was not found"
    
    # If found, we assume the molecule is an alpha-amino-acid zwitterion.
    return True, "Found the alpha-amino-acid zwitterion signature: an alpha carbon attached to [NH3+] and C(=O)[O-]"