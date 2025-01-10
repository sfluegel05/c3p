"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is characterized by complex, repetitive structural units, high molecular weight,
    and specific types of linkages.

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
    
    # Consider excluding entities with very low molecular weight first
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 500:  # Increased threshold to better capture polymers
        return False, "Molecular weight too low for a typical polymer"
    
    # Count heteroatoms (non-carbon and non-hydrogen) to check for complexity
    heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [6, 1])
    if heteroatom_count < 2:
        return False, "Too few heteroatoms to suggest polymer complexity"
    
    # Look for a high number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Insufficient rotatable bonds to be a polymer"

    # Look for various common functional group linkages that indicate repeated units
    repeating_unit_patterns = [
        Chem.MolFromSmarts("C=C"),
        Chem.MolFromSmarts("C-O-C"),  # Ethers, common in polyethers
        Chem.MolFromSmarts("O-C(=O)-C"),  # Esters, common in polyesters
        Chem.MolFromSmarts("N-C(=O)-C"),  # Amide linkages, common in polyamides
    ]
    
    for pattern in repeating_unit_patterns:
        repeating_matches = mol.GetSubstructMatches(pattern)
        if len(repeating_matches) >= 3:
            return True, f"Contains repeating units indicative of polymer structure: {pattern.GetSmarts()}"

    return False, "Does not match typical polymer characteristics"