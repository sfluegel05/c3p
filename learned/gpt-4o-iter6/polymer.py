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
    
    # Consider excluding entities with low molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 500:  # Reduced threshold for polymer-like behavior
        return False, "Molecular weight too low for a typical polymer-like entity"
    
    # Check for presence of significant repeating motifs
    linkages = ["C-C", "C=C", "C-O-C", "O=C-O", "N-C(=O)-C"]
    found_repeats = 0

    for pattern in linkages:
        if len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))) > 1:
            found_repeats += 1
    
    if found_repeats == 0:
        return False, "No significant repeating unit patterns found"
    
    # Check for a moderate number of rotatable bonds (inference of some flexibility)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Insufficient rotatable bonds to be considered a typical polymer"
    
    return True, "Contains characteristics of polymer structure including linkages and potential repetition"