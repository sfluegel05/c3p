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

    # Add hydrogens to correctly identify explicit hydrogens
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns for glycerol backbone
    # Glycerol backbone with three carbons and three hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)C(O)")
    matches = mol.GetSubstructMatches(glycerol_pattern)

    if not matches:
        return False, "No glycerol backbone found"

    # For each glycerol backbone found
    for match in matches:
        c_atoms = [mol.GetAtomWithIdx(idx) for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        o_atoms = [mol.GetAtomWithIdx(idx) for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]

        # Map each hydroxyl oxygen to its carbon
        hydroxyls = []
        for o_atom in o_atoms:
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor in c_atoms:
                    hydroxyls.append((neighbor, o_atom))
                    break

        # Keep track of substituents
        substituted_positions = 0
        substituent_types = []

        for c_atom, o_atom in hydroxyls:
            # Check if oxygen is bonded to any atom other than the carbon and hydrogen
            substituent_found = False
            for bond in o_atom.GetBonds():
                neighbor = bond.GetOtherAtom(o_atom)
                if neighbor.GetIdx() != c_atom.GetIdx():
                    if neighbor.GetAtomicNum() != 1:
                        substituent_found = True

                        # Determine substituent type
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            # Ether linkage
                            if neighbor.GetAtomicNum() == 6:
                                # Check for vinyl ether (alk-1-enyl)
                                is_vinyl = False
                                for next_bond in neighbor.GetBonds():
                                    nbr = next_bond.GetOtherAtom(neighbor)
                                    if nbr.GetIdx() != o_atom.GetIdx() and next_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                        is_vinyl = True
                                        break
                                if is_vinyl:
                                    substituent_types.append('alk-1-enyl')
                                else:
                                    substituent_types.append('alkyl')
                            else:
                                substituent_types.append('other')
                        elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            # Ester linkage (acyl)
                            substituent_types.append('acyl')
                        else:
                            substituent_types.append('other')
                    else:
                        # Oxygen bonded to hydrogen (hydroxyl group)
                        continue

            if substituent_found:
                substituted_positions += 1

        if substituted_positions == 2 and all(sub_type in ['acyl', 'alkyl', 'alk-1-enyl'] for sub_type in substituent_types):
            return True, f"Molecule is a diradylglycerol with {substituted_positions} substituted positions"
        else:
            continue  # Check next match if available

    return False, "Does not meet criteria for diradylglycerol"