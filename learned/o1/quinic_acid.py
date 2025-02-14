"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid or its derivative based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid - a cyclohexane ring with multiple hydroxyl groups (or their derivatives) and a carboxylic acid group attached to the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid or derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for cyclohexane ring
    cyclohexane_smarts = "[C&R]1[C&R][C&R][C&R][C&R][C&R]1"
    cyclohexane = Chem.MolFromSmarts(cyclohexane_smarts)
    ring_matches = mol.GetSubstructMatches(cyclohexane)
    if not ring_matches:
        return False, "No cyclohexane ring found"

    # For each cyclohexane ring found
    for match in ring_matches:
        ring_atoms = set(match)
        substituent_count = 0
        carboxylic_acid_found = False

        # Check substituents on ring atoms
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in ring_atoms:
                    # Check for carboxylic acid or ester group attached to ring
                    if neighbor.GetAtomicNum() == 6:
                        neighbor_atom = neighbor
                        # Look for -C(=O)O pattern
                        ester_carboxy_smarts = Chem.MolFromSmarts("C(=O)[O,N]")
                        bond = mol.GetBondBetweenAtoms(idx, neighbor_idx)
                        if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            frag = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
                            frag_smiles = Chem.MolToSmiles(frag, rootedAtAtom=neighbor_idx)
                            frag_mol = Chem.MolFromSmiles(frag_smiles)
                            if frag_mol is not None and frag_mol.HasSubstructMatch(ester_carboxy_smarts):
                                # Check if it's a carboxylic acid group
                                if neighbor_atom.GetDegree() == 3:
                                    carboxylic_acid_found = True
                                else:
                                    substituent_count += 1
                                continue

                    # Check for hydroxyl or ester group (-O-)
                    if neighbor.GetAtomicNum() == 8:
                        # Oxygen atom
                        # Check if it's a hydroxyl group
                        if neighbor.GetDegree() == 1:
                            substituent_count += 1
                            continue
                        elif neighbor.GetDegree() == 2:
                            # Possible ether or ester linkage
                            substituent_count +=1
                            continue

        if substituent_count >= 3 and carboxylic_acid_found:
            return True, "Contains cyclitol carboxylic acid core (quinic acid or derivative)"
        else:
            if substituent_count < 3:
                return False, f"Found {substituent_count} hydroxyl or ester groups on ring, need at least 3"
            if not carboxylic_acid_found:
                return False, "No carboxylic acid group attached to ring"

    return False, "Does not contain the required cyclitol carboxylic acid core"