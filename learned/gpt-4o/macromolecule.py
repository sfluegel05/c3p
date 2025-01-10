"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is characterized by a high molecular weight and the presence of repeating units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)
    # Setting a lower threshold for molecular weight
    if mol_wt < 1000:  # Reduced threshold to catch lighter macromolecules
        return False, f"Molecular weight too low for macromolecule: {mol_wt} Da"
    
    # Check for repeating units with more complex patterns
    # Using representative structures for biological macromolecules like peptides or sugars
    repeat_patterns = [
        Chem.MolFromSmarts("C(=O)N"),  # Peptide bond
        Chem.MolFromSmarts("C(=O)O"),  # Ester linkage
        Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C1"),  # Simple sugar ring, like glucose
    ]
    
    for idx, pattern in enumerate(repeat_patterns):
        matches = mol.GetSubstructMatches(pattern)
        if len(matches) >= 5:  # Requiring at least 5 matches to indicate repetition
            return True, f"Contains repeating units: Found {len(matches)} of pattern {idx}"
    
    return False, "No repeated units detected"

# Example usage
print(is_macromolecule("CC(=O)NCC(=O)NCC(=O)NCC(=O)O"))  # Example input