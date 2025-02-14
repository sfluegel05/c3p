"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: ribonucleoside
"""

from rdkit import Chem

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

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the ribose sugar ring (allowing modifications)
    # Ribose ring: five-membered ring with one oxygen and four carbons
    ribose_smarts = 'C1OC[C@H]([*])[C@@H]1[*]'
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)

    # Find matches for the ribose ring
    ribose_matches = mol.GetSubstructMatches(ribose_pattern, useChirality=False)
    if not ribose_matches:
        return False, "No ribose sugar ring found"

    # Define nucleobase patterns (purine and pyrimidine bases and their derivatives)
    nucleobase_patterns = [
        # Purine base pattern (adenine, guanine, and modifications)
        Chem.MolFromSmarts('n1c[nH]c2c1ncnc2'),
        Chem.MolFromSmarts('n1cnc2c1ncnc2'),
        Chem.MolFromSmarts('n1c[nH]c2c1nc[nH]c2'),
        # Pyrimidine base pattern (cytosine, uracil, thymine, and modifications)
        Chem.MolFromSmarts('c1c[nH]c(=O)[nH]c1=O'),
        Chem.MolFromSmarts('c1cnc(=O)[nH]c1=O'),
        Chem.MolFromSmarts('c1c[nH]c(=S)[nH]c1=O'),
        Chem.MolFromSmarts('c1c[nH]c(=O)[nH]c1=S'),
    ]

    # Check for nucleobase attachment via N-glycosidic bond
    nucleobase_found = False
    for match in ribose_matches:
        ribose_atoms = list(match)
        # Identify the anomeric carbon (connected to the oxygen in the ring)
        ring = mol.GetRingInfo().AtomRings()
        anomeric_carbons = []
        for ring_atoms in ring:
            if set(ribose_atoms).issubset(set(ring_atoms)):
                for idx in ring_atoms:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() == 'C':
                        neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
                        if 'O' in neighbors and 'N' in neighbors:
                            anomeric_carbons.append(idx)
                break  # We found the ribose ring

        if not anomeric_carbons:
            continue

        # Check if the nitrogen connected to the anomeric carbon is part of a nucleobase
        for anomeric_c_idx in anomeric_carbons:
            anomeric_c = mol.GetAtomWithIdx(anomeric_c_idx)
            for neighbor in anomeric_c.GetNeighbors():
                if neighbor.GetSymbol() == 'N':
                    # Create a combined molecule of the nucleobase
                    neighbor_idx = neighbor.GetIdx()
                    frag_atoms = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
                    base_indices = [i for i, frag in enumerate(frag_atoms) if neighbor_idx in frag]
                    base_atoms = frag_atoms[base_indices[0]]
                    base = Chem.PathToSubmol(mol, base_atoms)
                    # Check if the base matches any nucleobase pattern
                    for base_pattern in nucleobase_patterns:
                        if base.HasSubstructMatch(base_pattern):
                            nucleobase_found = True
                            break
                if nucleobase_found:
                    break
            if nucleobase_found:
                break
        if nucleobase_found:
            break

    if not nucleobase_found:
        return False, "No nucleobase attached via N-glycosidic bond"

    # Optional: Check for D-ribose configuration by looking at stereochemistry
    # Since stereochemistry might not always be specified, we relax this requirement

    return True, "Contains ribose sugar moiety attached to nucleobase via N-glycosidic bond"