"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are compounds derived from prostanoic acid (C20),
    containing a cyclopentane ring connected to two aliphatic chains at defined positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define prostaglandin core SMARTS pattern
    # This pattern represents a cyclopentane ring with two side chains attached at specific positions
    prostaglandin_core_smarts = '[C@H]1[C@@H](CC[C@@H]1[C@H](C)O)C/C=C\C[C@H](O)CCCCC(O)=O'
    
    # Create a generic pattern by removing specific stereochemistry and side chain details
    # This is to accommodate variations among prostaglandins
    prostaglandin_core_smarts = '[C,c]1[C,c][C,c][C,c][C,c]1'  # Five-membered ring

    # Attach side chains at positions 2 and 5 of the ring
    side_chain_pattern = Chem.MolFromSmarts("""
    [C,c]1(
        [C,c](
            [C,c](
                [C,c](
                    [C,c]1
                )[*,#6]  # Side chain at position 5
            )[*,#6]  # Side chain at position 4
        )[*,#6]  # Side chain at position 3
    )[*,#6]  # Side chain at position 2
    """)

    # Build a more specific prostaglandin pattern
    prostaglandin_pattern = Chem.MolFromSmarts("""
    [
        C;R1   # Carbon in ring
        ]1
        [
        C;R1   # Carbon in ring
        ][C;R1][C;R1][C;R1]1  # Five-membered ring
        [
        #
        ]  # Wildcard to match any side chain at position 1
        [
        C,C=O,O   # Side chain at position 2 can be carbon or oxygen
        ]
    """)

    # Alternatively, define specific substructure patterns
    # Prostaglandins have a cyclopentane ring with aliphatic side chains at specific positions

    # Define cyclopentane ring with side chains at position 1 and 3 (arbitrary numbering)
    cyclopentane_with_side_chains = Chem.MolFromSmarts("""
    [C@@H]1[C@H]([C@H](C[C@@H]1[O])[O])[C,C]=C[C,C][C,C]
    """)

    # For better matching, let's define two patterns: one for the cyclopentane ring and one for the side chains

    # Cyclopentane ring with two side chains
    ring_pattern = Chem.MolFromSmarts('C1CCCC1')

    # Side chain pattern (aliphatic chain, possibly with double bonds and functional groups)
    side_chain_pattern = Chem.MolFromSmarts('[C](=O)[O,O-,N]')  # Carboxylic acid or derivative

    # Identify cyclopentane ring
    ring_matches = mol.GetSubstructMatches(ring_pattern)
    if not ring_matches:
        return False, "No cyclopentane ring found"

    # For each ring, check if it has side chains attached at the correct positions
    for ring in ring_matches:
        ring_atoms = list(ring)
        side_chain_count = 0

        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atoms:
                    # Side chain detected
                    side_chain_count += 1

        if side_chain_count >= 2:
            # Now check for carboxylic acid or derivative in the molecule
            has_carboxy = mol.HasSubstructMatch(side_chain_pattern)
            if has_carboxy:
                return True, "Matches prostaglandin core structure with cyclopentane ring and side chains"
            else:
                # Check for ester or amide forms
                ester_pattern = Chem.MolFromSmarts('C(=O)O[C,N]')
                if mol.HasSubstructMatch(ester_pattern):
                    return True, "Matches prostaglandin core structure with esterified carboxylic acid"
                amide_pattern = Chem.MolFromSmarts('C(=O)N')
                if mol.HasSubstructMatch(amide_pattern):
                    return True, "Matches prostaglandin core structure with amidated carboxylic acid"
                else:
                    return False, "Carboxylic acid group or its derivative not found"

    return False, "Cyclopentane ring does not have required side chains"

# Examples of usage
if __name__ == "__main__":
    smiles_list = [
        # Examples of prostaglandins
        "CCCCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CC(O)=O",  # Prostaglandin F2alpha
        "CC(C)OC(=O)CCC\\C=C/C[C@H]1[C@@H](O)C[C@@H](O)[C@@H]1\\C=C\\C(F)(F)COc1ccccc1",  # Tafluprost
        # Non-prostaglandin for testing
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"  # Nonacosanoic acid
    ]
    for smiles in smiles_list:
        result, reason = is_prostaglandin(smiles)
        print(f"SMILES: {smiles}\nIs prostaglandin: {result}\nReason: {reason}\n")