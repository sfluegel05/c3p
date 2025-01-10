"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is characterized by repeating structural units, varied length, and common linkage groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of significant repeating motifs
    # This includes common polymer linkages
    linkage_patterns = [
        Chem.MolFromSmarts("C-C"),
        Chem.MolFromSmarts("C=C"),
        Chem.MolFromSmarts("C-O-C"),
        Chem.MolFromSmarts("O=C-O"),
        Chem.MolFromSmarts("N-C(=O)-C"),
        Chem.MolFromSmarts("C=C-C"),
        Chem.MolFromSmarts("C#C"),  # potentially capturing conjugated systems
    ]
    found_repeats = 0
    for pattern in linkage_patterns:
        if len(mol.GetSubstructMatches(pattern)) > 2:  # at least three instances of the pattern to imply repetition
            found_repeats += 1

    if found_repeats == 0:
        return False, "No significant repeating unit patterns found"
    
    # Check for moderate chain length
    # Use the number of atoms as an indicator of longer chain lengths typical in polymers
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 30:  # Reduced threshold to account for smaller polymers
        return False, "Insufficient molecular size for a characteristic polymer"

    # Check for a moderate to high number of rotatable bonds (inference of some flexibility)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Insufficient rotatable bonds to be considered a typical polymer"
    
    return True, "Contains characteristics of polymer structure including linkages and potential repetition"