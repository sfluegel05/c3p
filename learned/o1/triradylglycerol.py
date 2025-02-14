"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    glycerol_pattern = Chem.MolFromSmarts("""
    [C@@H](O*)-C(O*)-C(O*)
    """)
    if glycerol_pattern is None:
        return False, "Invalid glycerol backbone SMARTS pattern"

    # Search for the glycerol backbone
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No substituted glycerol backbone found"

    # Define patterns for acyl, alkyl, and alk-1-enyl groups
    acyl_pattern = Chem.MolFromSmarts("O=C-[#6]")
    alkyl_pattern = Chem.MolFromSmarts("[#6]")
    alk1enyl_pattern = Chem.MolFromSmarts("C=C-[#6]")

    # For each match, check the substituents
    for match in matches:
        is_valid = True
        reasons = []
        # Get the indices of the three carbons in the glycerol backbone
        c1_idx, c2_idx, c3_idx = match

        # List of carbons to check
        glycerol_carbons = [c1_idx, c2_idx, c3_idx]

        for idx in glycerol_carbons:
            carbon = mol.GetAtomWithIdx(idx)
            # Find the oxygen connected to this carbon
            oxygens = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetBondBetweenAtoms(nbr.GetIdx(), carbon.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE]
            if not oxygens:
                is_valid = False
                reasons.append(f"Carbon at index {idx} is not connected to an oxygen atom")
                continue
            oxygen = oxygens[0]
            # Find the substituent connected to the oxygen (excluding the glycerol carbon)
            substituents = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetIdx() != carbon.GetIdx()]
            if not substituents:
                is_valid = False
                reasons.append(f"Oxygen at index {oxygen.GetIdx()} has no substituent")
                continue
            substituent_atom = substituents[0]
            # Create a fragment of the substituent for pattern matching
            substituent_indices = Chem.rdmolops.GetShortestPath(mol, oxygen.GetIdx(), substituent_atom.GetIdx())
            substituent = Chem.PathToSubmol(mol, substituent_indices)
            # Check if substituent matches any of the patterns
            if substituent.HasSubstructMatch(acyl_pattern):
                reasons.append(f"Position {idx}: acyl substituent found")
            elif substituent.HasSubstructMatch(alk1enyl_pattern):
                reasons.append(f"Position {idx}: alk-1-enyl substituent found")
            elif substituent.HasSubstructMatch(alkyl_pattern):
                reasons.append(f"Position {idx}: alkyl substituent found")
            else:
                is_valid = False
                reasons.append(f"Position {idx}: substituent is not acyl, alkyl, or alk-1-enyl")
        if is_valid:
            return True, "Molecule is a triradylglycerol"
        else:
            return False, "; ".join(reasons)

    return False, "No valid triradylglycerol structure found"