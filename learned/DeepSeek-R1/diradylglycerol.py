"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:64435 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol has a glycerol backbone with exactly two substituent groups
    (acyl, alkyl, or alk-1-enyl) attached via oxygen atoms.

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

    # Glycerol backbone pattern: three contiguous carbons (C-C-C)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Check each possible glycerol backbone configuration
    for backbone in matches:
        substituent_count = 0
        # Iterate through each carbon in the glycerol backbone
        for carbon_idx in backbone:
            carbon = mol.GetAtomWithIdx(carbon_idx)
            # Check all oxygen atoms attached to this carbon
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() != 8:
                    continue  # Not an oxygen atom
                oxygen_idx = neighbor.GetIdx()
                # Ensure the bond is a single bond (ester/ether/vinyl ether)
                bond = mol.GetBondBetweenAtoms(carbon_idx, oxygen_idx)
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Get the atom attached to the oxygen (substituent)
                substituent_atoms = []
                for o_neighbor in neighbor.GetNeighbors():
                    if o_neighbor.GetIdx() != carbon_idx:
                        substituent_atoms.append(o_neighbor)
                # Check each substituent atom (typically one)
                for sub_atom in substituent_atoms:
                    # Acyl (ester): substituent is C connected to O via C=O
                    if sub_atom.GetAtomicNum() == 6:
                        # Check if this carbon is part of a carbonyl group
                        for bond in sub_atom.GetBonds():
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                other_atom = bond.GetOtherAtom(sub_atom)
                                if other_atom.GetAtomicNum() == 8:
                                    substituent_count += 1
                                    break
                        else:
                            # Alkyl (ether): check for at least one carbon chain
                            # Alk-1-enyl (vinyl ether): check for adjacent double bond
                            has_double = any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in sub_atom.GetBonds())
                            if has_double:
                                # Check if the double bond is between C1 and C2 (alk-1-enyl)
                                next_atoms = [a for a in sub_atom.GetNeighbors() if a.GetIdx() != oxygen_idx]
                                for a in next_atoms:
                                    if mol.GetBondBetweenAtoms(sub_atom.GetIdx(), a.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                        substituent_count += 1
                                        break
                            else:
                                # Assume alkyl if no double bonds and has at least one more carbon
                                if sub_atom.GetDegree() >= 1:
                                    substituent_count += 1
                    # If substituent is hydrogen (hydroxyl group), ignore
        # Check if exactly two substituents found
        if substituent_count == 2:
            return True, "Contains glycerol backbone with exactly two substituent groups (acyl/alkyl/alk-1-enyl)"
    
    return False, "Found fewer or more than two valid substituent groups"