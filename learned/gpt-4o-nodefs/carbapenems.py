"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem typically features a bicyclic ring structure with a β-lactam
    fused to a five-membered ring, often containing a sulfur atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a refined SMARTS pattern for carbapenem
    # β-lactam ring with nitrogen, carbonyl group, fused to five-membered ring with sulfur
    # Carbons in β-lactam are typically sp3, and stereochemistry may be specified
    carbapenem_pattern = Chem.MolFromSmarts("C1N[C@@H]2[C@@H](C)C(=O)N2C1O")  # Example pattern, check with specifics

    if not mol.HasSubstructMatch(carbapenem_pattern):
        return False, "Does not match the carbapenem structural motif"

    # Check for sulfur in the structure, often found in many carbapenems
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if sulfur_count == 0:
        return False, "No sulfur atom found, which is less likely for a typical carbapenem"

    return True, "Contains the β-lactam core with a fused ring and potential sulfur, typical of carbapenems"

# Further refinements can be made based on a more precise understanding of the defining structural features of carbapenems