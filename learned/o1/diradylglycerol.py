"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol molecule where two of the three hydroxyl groups
    are substituted with acyl (ester-linked), alkyl (ether-linked), or alk-1-enyl (vinyl ether-linked) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")  # Simplified glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Find all glycerol backbones in the molecule
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Check each glycerol backbone in the molecule
    for match in matches:
        # Get the atoms corresponding to the glycerol carbons
        glycerol_atoms = list(match)
        # Initialize counts for substituted and free hydroxyl groups
        substituted_positions = 0
        free_hydroxyl_positions = 0

        # For each position (carbon atom in glycerol)
        for idx in glycerol_atoms:
            atom = mol.GetAtomWithIdx(idx)
            # Get all bonds from this atom
            bonds = atom.GetBonds()
            # Check for connections that indicate substitution
            is_substituted = False
            for bond in bonds:
                neighbor = bond.GetOtherAtom(atom)
                # Skip hydrogens
                if neighbor.GetAtomicNum() == 1:
                    continue
                # Check if the neighbor is oxygen
                if neighbor.GetAtomicNum() == 8:
                    # Check the bond from oxygen to its substituent
                    oxygen_atom = neighbor
                    oxygen_bonds = oxygen_atom.GetBonds()
                    # If oxygen is connected to something other than the glycerol carbon and hydrogen
                    for o_bond in oxygen_bonds:
                        o_neighbor = o_bond.GetOtherAtom(oxygen_atom)
                        if o_neighbor.GetIdx() != atom.GetIdx() and o_neighbor.GetAtomicNum() != 1:
                            # Check for ester or ether linkage
                            if o_neighbor.GetAtomicNum() == 6:
                                # Check if it's part of ester or ether
                                if o_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                    is_substituted = True
                                    break
                else:
                    # Direct substitution (unlikely in glycerol backbone)
                    is_substituted = True
            if is_substituted:
                substituted_positions += 1
            else:
                free_hydroxyl_positions += 1

        # Check if exactly two positions are substituted
        if substituted_positions == 2:
            return True, "Molecule is a diradylglycerol with two substituted positions"
        else:
            return False, f"Found {substituted_positions} substituted positions, expected 2"

    return False, "No matching glycerol backbone with two substitutions found"