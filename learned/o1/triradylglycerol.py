"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol
"""

from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define the SMARTS pattern for the glycerol backbone with substitutions
    # The pattern represents a glycerol backbone where each carbon is connected to an oxygen,
    # and each oxygen is attached to a substituent (any group except hydrogen)
    glycerol_pattern = Chem.MolFromSmarts("""
    [C@@H]([O][!H0])[C@@H]([O][!H0])[C@@H]([O][!H0])
    """)
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No substituted glycerol backbone found"

    # Define SMARTS patterns for acyl, alkyl, and alk-1-enyl groups
    acyl_pattern = Chem.MolFromSmarts("C(=O)[C,H]")
    alkyl_pattern = Chem.MolFromSmarts("[C][C,H]")
    alk1enyl_pattern = Chem.MolFromSmarts("C=C[C,H]")

    # Iterate through possible glycerol backbones
    for match in matches:
        is_triradyl = True
        reasons = []
        for idx in match:
            oxygen = None
            atom = mol.GetAtomWithIdx(idx)
            # Find the oxygen connected to this carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    oxygen = neighbor
                    break
            if oxygen is None:
                is_triradyl = False
                reasons.append(f"Carbon at index {idx} has no connected oxygen")
                break

            # Check the substituent connected to the oxygen
            substituents = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetIdx() != atom.GetIdx()]
            if len(substituents) != 1:
                is_triradyl = False
                reasons.append(f"Oxygen at index {oxygen.GetIdx()} does not have exactly one substituent")
                break

            substituent_atom = substituents[0]
            # Create a fragment of the substituent for substructure matching
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, substituent_atom.GetIdx())
            submol = Chem.PathToSubmol(mol, env)
            # Check if substituent is acyl, alkyl, or alk-1-enyl
            if submol.HasSubstructMatch(acyl_pattern):
                reasons.append(f"Position {idx}: acyl substituent found")
            elif submol.HasSubstructMatch(alk1enyl_pattern):
                reasons.append(f"Position {idx}: alk-1-enyl substituent found")
            elif submol.HasSubstructMatch(alkyl_pattern):
                reasons.append(f"Position {idx}: alkyl substituent found")
            else:
                is_triradyl = False
                reasons.append(f"Position {idx}: substituent is not acyl, alkyl, or alk-1-enyl")
                break

        if is_triradyl:
            return True, "Molecule is a triradylglycerol"
        else:
            continue  # Try next possible glycerol backbone

    return False, "; ".join(reasons) if reasons else "No valid triradylglycerol structure found"