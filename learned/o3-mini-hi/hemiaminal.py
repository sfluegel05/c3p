"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemi-aminal compounds
Definition: Any organic amino compound that has an amino group and a hydroxy group 
attached to the same carbon atom. Hemiaminals are intermediates in the formation of imines.
Examples include alpha-hydroxyglycine, isopseudostrychnine, and others in the provided list.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is defined as a molecule containing a carbon atom bearing both a 
    hydroxyl group (–OH) and an amino group (–NH, –NHR, or –NR2) attached to the same carbon.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule contains a hemiaminal motif, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the hemiaminal group.
    # This pattern looks for a tetrahedral (sp3) carbon that is not part of a carbonyl (C=O)
    # and that is attached to both an -OH group ([OX2H]) and an amine group ([NX3]).
    hemiaminal_pattern = Chem.MolFromSmarts("[C;!$(C=O)]([OX2H])([NX3])")
    
    # Check if the molecule has at least one match for the hemiaminal pattern.
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Molecule contains a hemiaminal motif (a carbon atom with both -OH and -N groups)."
    else:
        return False, "No hemiaminal motif (a carbon atom bearing both -OH and -N groups) was found."

# Example usage (uncomment the lines below to test):
# test_smiles = "NC(O)C(O)=O"  # Example: alpha-hydroxyglycine
# result, reason = is_hemiaminal(test_smiles)
# print(result, reason)