"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:18385 alditol
"""

from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH,
    formally derivable from an aldose by reduction of the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for an alditol chain
    # This matches an acyclic chain of carbons with attached hydroxyl groups,
    # starting and ending with CH2OH groups.
    pattern_str = '[O][CH2][CH](O){1,10}[CH2][O]'
    pattern = Chem.MolFromSmarts(pattern_str)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for matches in the molecule
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "Molecule does not contain an alditol chain"

    # Check each match for additional criteria
    for match in matches:
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]

        # Check that the matched atoms are not in rings
        if any(atom.IsInRing() for atom in atoms):
            continue  # Skip this match

        # Check that all carbons are sp3 hybridized and have OH groups
        carbons = [atom for atom in atoms if atom.GetAtomicNum() == 6]
        all_sp3 = all(atom.GetHybridization() == Chem.HybridizationType.SP3 for atom in carbons)
        all_has_OH = True
        for atom in carbons:
            # Check that carbon has an attached hydroxyl group
            has_OH = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        has_OH = True
                        break
            if not has_OH:
                all_has_OH = False
                break
        if not all_sp3 or not all_has_OH:
            continue  # Skip this match

        # Passed all checks
        n_carbons = len(carbons)
        return True, f"Contains an alditol chain of length {n_carbons} carbons"

    # No valid alditol chain found
    return False, "Molecule does not contain an alditol chain"