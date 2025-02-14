"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol molecule bearing a single acyl, alkyl, or alk-1-enyl substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize molecule
    Chem.SanitizeMol(mol)

    # Define glycerol pattern (glycerol backbone with three carbons and three hydroxyl groups)
    glycerol_smarts = "[OX2H][CX4H]([OX2H])[CX4H2][OX2H]"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    matches = mol.GetSubstructMatches(glycerol_pattern)

    if not matches:
        return False, "No glycerol backbone found"

    # Assume the first match is the glycerol backbone
    glycerol_match = matches[0]
    glycerol_atoms = set(glycerol_match)

    # Count the number of substituents on the glycerol carbons
    substituent_count = 0
    substituent_types = []

    for atom_idx in glycerol_match:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            continue  # Skip non-carbon atoms
        neighbors = atom.GetNeighbors()
        for neighbor in neighbors:
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in glycerol_atoms:
                substituent_count += 1
                # Check substituent type
                substituent_atom = neighbor

                # Acyl group: O=C-R (ester linkage)
                acyl_smarts = "[CX3](=O)[O][CX4]"
                acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
                acyl_match = mol.GetSubstructMatches(acyl_pattern)
                if acyl_match:
                    substituent_types.append("acyl")
                    continue

                # Alkyl group: simple carbon chain
                # Check if substituent is an alkyl chain
                alkyl_smarts = "[CX4H2][CX4H2]"
                alkyl_pattern = Chem.MolFromSmarts(alkyl_smarts)
                alkyl_match = mol.GetSubstructMatches(alkyl_pattern)
                if alkyl_match:
                    substituent_types.append("alkyl")
                    continue

                # Alk-1-enyl group: C=C-C
                alk1enyl_smarts = "[CX3]=[CX3][CX4]"
                alk1enyl_pattern = Chem.MolFromSmarts(alk1enyl_smarts)
                alk1enyl_match = mol.GetSubstructMatches(alk1enyl_pattern)
                if alk1enyl_match:
                    substituent_types.append("alk-1-enyl")
                    continue

                substituent_types.append("other")

    if substituent_count != 1:
        return False, f"Glycerol has {substituent_count} substituents, need exactly 1"

    substituent_type = substituent_types[0] if substituent_types else "unknown"

    if substituent_type not in ["acyl", "alkyl", "alk-1-enyl"]:
        return False, "Substituent is not acyl, alkyl, or alk-1-enyl group"

    return True, f"Molecule is a monoradylglycerol with one {substituent_type} substituent"