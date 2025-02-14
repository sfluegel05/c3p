"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    Focuses on identifying complex chains and repeating units, common in polymers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is identified as a polymer-like structure, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define polymer-like repeating unit patterns
    patterns = [
        Chem.MolFromSmarts("[C](=[O,N])-[C]"),  # polymethacrylate motifs
        Chem.MolFromSmarts("[C]([C])[C]([C])[C]"),  # generally for polyolefins
        Chem.MolFromSmarts("[O]-[C]-[O]"),  # linkage common in polyethers
        Chem.MolFromSmarts("C=C-C"),  # unsaturated backbone common in rubbers
        Chem.MolFromSmarts("C-C-C-C"),  # flexible hydrocarbon chain
    ]

    # Check matches for polymer-like chains
    total_matches = 0
    for pattern in patterns:
        matches = mol.GetSubstructMatches(pattern)
        total_matches += len(matches)

    # Lower threshold to identify polymers more consistently
    if total_matches > 5:  # refined threshold
        return True, f"Contains repeating substructural patterns: found {total_matches}"

    # Consider large size and patterns together for polymer detection
    if mol.GetNumAtoms() > 50 and total_matches > 2:  # balance of size and pattern detection
        return True, f"Large molecular size with repeating units, possible polymer ({mol.GetNumAtoms()} atoms)"

    # Check for long single bond chains (indicative of polymer backbone)
    dm = rdmolops.GetDistanceMatrix(mol)
    if dm.size > 0 and max(dm.flatten()) > 30:  # refined chain length threshold
        return True, "Contains long chains indicative of polymeric structure"

    return False, "No explicit polymer-like features detected"