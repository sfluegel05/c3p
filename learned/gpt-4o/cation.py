"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is defined as a monoatomic or polyatomic species having one or more elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate net charge of the molecule
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

    # Check if net charge is positive
    if net_charge > 0:
        # Additional checking for specific cation groups
        # Here we define known cation groups' SMARTS patterns
        cation_patterns = [
            '[N+]',     # Generic quaternary ammonium
            '[P+]',     # Phosphonium
            '[S+]',     # Sulfonium
            '[O+]',     # Oxonium (less common but valid)
            '[C+]',     # Carbocation
            '[Fe+3]',   # Ferric ion example
            # Add more if necessary
        ]
        
        for pattern in cation_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                return True, f"Positive charge of {pattern} detected"

    return False, "No cationic groups found"