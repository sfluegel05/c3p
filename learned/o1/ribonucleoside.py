"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: ribonucleoside
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the ribose sugar (allowing for common modifications)
    # Ribose ring pattern (five-membered ring with oxygen and four carbons)
    ribose_smarts = '[C@H]1([O])[C@@H]([O])[C@H]([O])[C@@H]1[O]'
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    if ribose_pattern is None:
        return False, "Error in ribose pattern"

    # Find matches for the ribose sugar
    ribose_matches = mol.GetSubstructMatches(ribose_pattern, useChirality=True)
    if not ribose_matches:
        return False, "No D-ribose sugar moiety found"

    # Define nucleobase patterns (allowing for common modifications)
    nucleobase_patterns = []

    # Purine base pattern (adenine, guanine, and modified purines)
    purine_smarts = 'n1c([nH])cnc1n'
    purine_pattern = Chem.MolFromSmarts(purine_smarts)
    nucleobase_patterns.append(purine_pattern)

    # Pyrimidine base pattern (cytosine, uracil, thymine, and modified pyrimidines)
    pyrimidine_smarts = 'c1cnc([nH])c1=O'
    pyrimidine_pattern = Chem.MolFromSmarts(pyrimidine_smarts)
    nucleobase_patterns.append(pyrimidine_pattern)

    # Check for nucleobase attached to sugar via N-glycosidic bond
    nucleobase_found = False
    for ribose_match in ribose_matches:
        ribose_atoms = set(ribose_match)
        anomeric_c = None
        for idx in ribose_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetDegree() == 3 and atom.GetSymbol() == 'C':
                anomeric_c = atom
                break
        if anomeric_c is None:
            continue

        # Check connections from anomeric carbon to nucleobase nitrogen
        connected = False
        for nbr in anomeric_c.GetNeighbors():
            if nbr.GetSymbol() == 'N':
                # Check if the nitrogen is part of a nucleobase
                for pattern in nucleobase_patterns:
                    if mol.HasSubstructMatch(pattern):
                        nucleobase_found = True
                        connected = True
                        break
            if connected:
                break
        if nucleobase_found:
            break

    if not nucleobase_found:
        return False, "No nucleobase moiety attached to D-ribose sugar via N-glycosidic bond"

    # Verify stereochemistry of D-ribose
    stereocenters = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    d_ribose_centers = [
        (idx, 'R') for idx, config in stereocenters if idx in ribose_atoms and config == 'R'
    ]
    if len(d_ribose_centers) < 3:
        return False, "Insufficient stereochemistry to confirm D-ribose"

    return True, "Contains D-ribose sugar moiety attached to nucleobase via N-glycosidic bond"