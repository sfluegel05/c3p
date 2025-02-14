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

    # Define a more specific glycerol backbone pattern with oxygen atoms
    glycerol_pattern = Chem.MolFromSmarts("""
    [C:1]-[C:2]-[C:3],
    [C:1][O],
    [C:2][O],
    [C:3][O]
    """)

    # Alternative SMARTS pattern for glycerol backbone with attached oxygens
    glycerol_pattern = Chem.MolFromSmarts("""
    [$([C@@H](O)-[CH2]-O),$([C@@H](O)-[CH2]-O)]
    """)

    # For this example, we'll use a pattern that matches glycerol backbone with oxygens
    glycerol_pattern = Chem.MolFromSmarts("[C;H2][C;H][C;H2]")

    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # For each glycerol backbone found, check for ether linkages
    for match in matches:
        c1, c2, c3 = [mol.GetAtomWithIdx(idx) for idx in match]

        # Check for oxygen atoms connected to each carbon
        c_atoms = [c1, c2, c3]
        ether_found = False

        for c in c_atoms:
            # Get all oxygen neighbors
            oxygen_neighbors = [nbr for nbr in c.GetNeighbors() if nbr.GetAtomicNum() == 8]

            for o in oxygen_neighbors:
                bond = mol.GetBondBetweenAtoms(c.GetIdx(), o.GetIdx())
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    continue  # Skip non-single bonds

                # Check if oxygen is connected to another carbon (not part of glycerol backbone)
                other_carbons = [nbr for nbr in o.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != c.GetIdx()]
                for oc in other_carbons:
                    # Ensure the other carbon is not a carbonyl carbon
                    is_carbonyl = False
                    for obond in oc.GetBonds():
                        if obond.GetBondType() == Chem.BondType.DOUBLE and obond.GetOtherAtom(oc).GetAtomicNum() == 8:
                            is_carbonyl = True
                            break
                    if not is_carbonyl:
                        # Found an ether linkage
                        ether_found = True
                        break
                if ether_found:
                    break
            if ether_found:
                break

        if ether_found:
            return True, "Contains glycerol backbone with at least one alkyl chain attached via ether linkage"

    return False, "No ether linkage found on glycerol backbone"