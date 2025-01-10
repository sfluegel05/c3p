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
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)
    # Adjusting the threshold to consider larger macromolecules
    if mol_wt < 2000:  # Increasing the threshold for molecular weight
        return False, f"Molecular weight too low for macromolecule: {mol_wt} Da"
    
    # Check for repeating units with more complex patterns
    # Using more representative structures for biological macromolecules like peptides or sugars
    repeat_patterns = [
        Chem.MolFromSmarts("C(=O)N"),  # Peptide bond
        Chem.MolFromSmarts("C(=O)OC"),  # Ester linkage often found in polysaccharides
        Chem.MolFromSmarts("CCO[C@H](O)[C@H](O)"),  # Example of sugar pattern
    ]
    
    for idx, pattern in enumerate(repeat_patterns):
        matches = mol.GetSubstructMatches(pattern)
        if len(matches) > 10:  # Requiring more than 10 matches to indicate repetition
            return True, f"Contains repeating units: Found {len(matches)} of pattern {idx}"
    
    return False, "No repeated units detected"

# Example usage:
print(is_macromolecule("CC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)N"))