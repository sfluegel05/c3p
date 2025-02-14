"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:64675 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is a lipid similar in structure to a glycerolipid but in which
    one or more of the carbon atoms on glycerol is bonded to an alkyl chain via an
    ether linkage, as opposed to the usual ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern (3 consecutive carbons with hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("C[C@H](O)CO")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Flag to check if at least one ether linkage is present
    has_ether_linkage = False

    # Iterate over glycerol backbone matches
    for match in matches:
        glycerol_carbons = match[:3]  # Indices of the three carbons in glycerol

        # Check each carbon in glycerol
        for idx in glycerol_carbons:
            atom = mol.GetAtomWithIdx(idx)
            oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]

            # Check bonds between carbon and oxygen
            for oxygen in oxygen_neighbors:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx())
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    # Check if oxygen is connected to an alkyl chain
                    oxygen_bonds = oxygen.GetBonds()
                    for obond in oxygen_bonds:
                        other_atom = obond.GetOtherAtom(oxygen)
                        if other_atom.GetIdx() != atom.GetIdx() and other_atom.GetAtomicNum() == 6:
                            # Found an ether linkage
                            has_ether_linkage = True
                            break

            # Break early if we have found an ether linkage
            if has_ether_linkage:
                break
        if has_ether_linkage:
            break

    if not has_ether_linkage:
        return False, "No ether linkage found on glycerol backbone"

    return True, "Contains glycerol backbone with at least one alkyl chain attached via ether linkage"