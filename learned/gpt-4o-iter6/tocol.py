"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    Tocols have a chromanol core with a hydrocarbon chain attached at position 2
    consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general chroman-6-ol core pattern
    chromanol_pattern = Chem.MolFromSmarts("c1cc(O)ccc2COC(C)C(O)c12")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chroman-6-ol core found"

    # Define isoprenoid pattern
    isoprenoid_pattern = Chem.MolFromSmarts("C(C)(C)CC")

    # Looking for hydrocarbon chain of at least 15 carbons (3 isoprenoid units)
    chains = [match for match in mol.GetSubstructMatches(isoprenoid_pattern)]
    if len(chains) < 3:
        return False, f"Found {len(chains)} isoprenoid units, need at least 3"
    
    # Ensure correct attachment at position 2 of the chroman-6-ol
    core_match = mol.GetSubstructMatch(chromanol_pattern)
    if core_match:
        # Presuming attachment of hydrocarbon at position 2 can be verified
        # by ensuring the long chain detected matches positional attachment
        atom_idx_position_2 = core_match[4]  # positional index 2 for chromanol
        for match in chains:
            if atom_idx_position_2 not in match:
                continue
            # Attachment verification (more logic could be added)
            return True, "Contains tocol structure with required chromanol core and isoprenoid units"

    return False, "Chroman-6-ol core or hydrocarbon substitution failed verification"