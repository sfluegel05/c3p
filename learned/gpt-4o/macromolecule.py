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
    if mol_wt < 1000:  # Arbitrarily choosing 1000 Da as a lower limit for macromolecules
        return False, f"Molecular weight too low for macromolecule: {mol_wt} Da"
    
    # Check for repeating units - this might be complex and is a heuristic
    # Looking for repeating ethylene or peptide-like substructures as an example
    repeat_patterns = [
        Chem.MolFromSmarts("C=C"),  # Ethylene
        Chem.MolFromSmarts("C(=O)N"),  # Peptide bond
    ]
    
    for pattern in repeat_patterns:
        matches = mol.GetSubstructMatches(pattern)
        if len(matches) > 5:  # Arbitrarily considering more than 5 matches as indicating repetition
            return True, f"Contains repeating units: Found {len(matches)} of pattern {pattern}"
    
    return False, "No repeated units detected"

# Example usage:
print(is_macromolecule("CC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)N"))