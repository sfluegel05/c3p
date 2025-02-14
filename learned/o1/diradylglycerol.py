"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem

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

    # Glycerol backbone with optional substitution at any position
    glycerol_pattern = Chem.MolFromSmarts("C([O])(C([O])C([O]))")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # For each glycerol backbone found
    for match in matches:
        # Indices of carbons and oxygens in the glycerol backbone
        c1_idx, o1_idx, c2_idx, o2_idx, c3_idx, o3_idx = match

        # Keep track of the number of substituted positions
        substituted_positions = 0

        # List to store types of substituents
        substituent_types = []

        # Function to check for substituents on each hydroxyl group
        def check_substitution(carbon_idx, oxygen_idx):
            carbon = mol.GetAtomWithIdx(carbon_idx)
            oxygen = mol.GetAtomWithIdx(oxygen_idx)

            # Bonds from oxygen
            for bond in oxygen.GetBonds():
                neighbor = bond.GetOtherAtom(oxygen)
                if neighbor.GetIdx() != carbon_idx:
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        # Check for ether linkage (alkyl or alk-1-enyl)
                        substituent = neighbor.GetIdx()
                        substituent_atom = mol.GetAtomWithIdx(substituent)
                        if substituent_atom.GetAtomicNum() == 6:
                            # Check for vinyl ether (alk-1-enyl)
                            if substituent_atom.GetDegree() > 1:
                                for nb in substituent_atom.GetNeighbors():
                                    if nb.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(substituent_atom.GetIdx(), nb.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                        substituent_types.append('alk-1-enyl')
                                        return 1
                            substituent_types.append('alkyl')
                            return 1
                    elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        # Check for ester linkage (acyl)
                        substituent_types.append('acyl')
                        return 1
            return 0

        # Check substitutions at each position
        substituted_positions += check_substitution(c1_idx, o1_idx)
        substituted_positions += check_substitution(c2_idx, o2_idx)
        substituted_positions += check_substitution(c3_idx, o3_idx)

        if substituted_positions >= 2:
            # Verify that substituents are acyl, alkyl, or alk-1-enyl
            if all(sub_type in ['acyl', 'alkyl', 'alk-1-enyl'] for sub_type in substituent_types):
                return True, f"Molecule is a diradylglycerol with {substituted_positions} substituted positions"
            else:
                return False, "Substituents are not acyl, alkyl, or alk-1-enyl groups"

    return False, "Does not meet criteria for diradylglycerol"