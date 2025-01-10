"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: CHEBI:36235 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid core conjugated to a functional group that adds hydrophilicity or charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general bile acid core pattern (steroid-like structure)
    # This pattern captures the essential 4-ring steroid structure with potential hydroxyl groups
    bile_acid_core_pattern = Chem.MolFromSmarts("[C@H]1[C@@]2([C@H]([C@@H]3[C@]([C@@H]([C@H]4[C@](CC3)(CC[C@H]4[OH0,OH1])[H])C)[H])CC[C@@]2([C@@H](C1)[OH0,OH1])[H])")
    if not mol.HasSubstructMatch(bile_acid_core_pattern):
        return False, "No bile acid core structure found"

    # Define patterns for conjugated groups
    conjugated_groups = [
        Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4H2]"),  # Glycine
        Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])[CX4H2][CX4H2][NX3]"),  # Taurine
        Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])([OX2])"),  # Sulfate
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])"),  # Glucuronic acid
        Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX3](=[OX1])"),  # Coenzyme A-like structure
        Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4H2][CX4H2][SX4](=[OX1])(=[OX1])"),  # Taurine variant
        Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX4H2][CX4H2][CX4H2][NX3]"),  # Amino acid variants
    ]

    # Check for presence of conjugated groups attached to the bile acid core
    has_conjugated_group = False
    for pattern in conjugated_groups:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # Check if the conjugated group is attached to the bile acid core
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in match:
                        # Check if neighbor is part of the bile acid core
                        if mol.HasSubstructMatch(bile_acid_core_pattern, useChirality=True, atomMap={neighbor.GetIdx(): 0}):
                            has_conjugated_group = True
                            break
                if has_conjugated_group:
                    break
            if has_conjugated_group:
                break
        if has_conjugated_group:
            break

    if not has_conjugated_group:
        return False, "No conjugated group found attached to bile acid core"

    # Check for hydrophilicity or charge
    # Calculate the number of polar atoms (O, N, S) and formal charges
    polar_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8, 16])
    formal_charges = sum(abs(atom.GetFormalCharge()) for atom in mol.GetAtoms())

    if polar_atoms < 3 and formal_charges == 0:
        return False, "Insufficient hydrophilicity or charge"

    return True, "Contains bile acid core with conjugated hydrophilic or charged group"