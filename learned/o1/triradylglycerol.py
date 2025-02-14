"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol
"""

from rdkit import Chem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol is a glycerol compound where each of the three hydroxyl groups 
    is substituted with either an acyl, alkyl, or alk-1-enyl group at positions 
    sn-1, sn-2, and sn-3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the substituted glycerol backbone
    # Pattern: glycerol backbone with three positions substituted via ester or ether linkages
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO[!H0])[C@@H](CO[!H0])[C@@H](CO[!H0])")
    if glycerol_pattern is None:
        return False, "Invalid glycerol backbone SMARTS pattern"

    # Find matches to the glycerol backbone pattern
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No substituted glycerol backbone found"

    # Define SMARTS patterns for acyl, alkyl, and alk-1-enyl groups attached via oxygen
    # Acyl group attached via ester linkage: -O-C(=O)C
    acyl_pattern = Chem.MolFromSmarts("COC(=O)[C;!H0]")
    # Alkyl group attached via ether linkage: -O-C-C
    alkyl_pattern = Chem.MolFromSmarts("CO[C;!H0][C;!H0]")
    # Alk-1-enyl group attached via ether linkage: -O-C=C-C
    alk1enyl_pattern = Chem.MolFromSmarts("COC=C[C;!H0]")

    if None in [acyl_pattern, alkyl_pattern, alk1enyl_pattern]:
        return False, "Invalid substituent SMARTS pattern"

    # Check each substituent at the glycerol backbone positions
    is_triradyl = True
    reasons = []
    for match in matches:
        is_triradyl = True
        reasons = []
        # Get the indices for the three carbon atoms in the glycerol backbone
        c1_idx, c2_idx, c3_idx = match[0], match[1], match[2]

        # List to store the oxygen atoms connected to each carbon
        oxygens = []
        for c_idx in [c1_idx, c2_idx, c3_idx]:
            carbon = mol.GetAtomWithIdx(c_idx)
            # Find oxygen atom connected to this carbon
            oxygen = None
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    oxygen = neighbor
                    break
            if oxygen is None:
                is_triradyl = False
                reasons.append(f"Carbon at index {c_idx} has no connected oxygen")
                break
            oxygens.append(oxygen)

        if not is_triradyl:
            continue  # Move to next match if any carbon lacks an oxygen

        # Check substituents attached to each oxygen
        substituent_patterns = [acyl_pattern, alkyl_pattern, alk1enyl_pattern]
        for oxygen in oxygens:
            # Get substituents attached to the oxygen (excluding the glycerol carbon)
            substituents = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetAtomicNum() != 6 or nbr.GetIdx() not in match]
            if len(substituents) != 1:
                is_triradyl = False
                reasons.append(f"Oxygen at index {oxygen.GetIdx()} does not have exactly one substituent")
                break

            # Create a fragment for the substituent
            substituent_atom = substituents[0]
            paths = Chem.rdmolops.GetShortestPath(mol, oxygen.GetIdx(), substituent_atom.GetIdx())
            substituent_frag = Chem.PathToSubmol(mol, paths)

            # Check if substituent matches any of the allowed patterns
            if any(substituent_frag.HasSubstructMatch(pattern) for pattern in substituent_patterns):
                continue
            else:
                is_triradyl = False
                reasons.append(f"Substituent at oxygen index {oxygen.GetIdx()} is not acyl, alkyl, or alk-1-enyl")
                break

        if is_triradyl:
            return True, "Molecule is a triradylglycerol"
        else:
            continue  # Try the next match

    return False, "; ".join(reasons) if reasons else "No valid triradylglycerol structure found"