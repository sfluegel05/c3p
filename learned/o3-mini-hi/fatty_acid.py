"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: fatty acid
Definition: “Any aliphatic monocarboxylic acid.
Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched and even‐numbered),
which may be saturated or unsaturated. By extension, the term is sometimes used to embrace all acyclic aliphatic carboxylic acids.”
This classifier now only requires that a free carboxyl (–C(=O)OH) exists and that at least one candidate contiguous carbon chain
attached to its carbonyl carbon is unbranched, does not include additional carbonyl (C=O) bonds or ring atoms, 
and has a total length (including the acid carbon) between 4 and 28.
Note: We no longer reject molecules just because they contain rings outside the candidate chain.
Also, we disallow chains that have C=O (other than the acid’s carbonyl) as these indicate extra functionalities.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if the given SMILES is a fatty acid based on:
      - Being a valid molecule,
      - Containing a free (undissociated) carboxylic acid group –C(=O)[OH],
      - Having at least one candidate contiguous, unbranched carbon chain (starting at the carbonyl carbon) that,
        when including the acid carbon, is between 4 and 28 carbons in length,
      - And the candidate chain itself does not include extra carbonyl functionalities or ring atoms.

    Args:
      smiles (str): SMILES string for the molecule.

    Returns:
      bool: True if the molecule is classified as a fatty acid, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Rule: Do not allow nitrogen atoms which signal peptides/amino acid derivatives.
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Molecule contains nitrogen atoms indicating peptide or amino acid derivatives"
    
    # Rule: Look for free carboxyl (–C(=O)[OH]) group.
    acid_smarts = "[CX3](=O)[OX2H]"  # free acid signature
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No free carboxyl (–C(=O)[OH]) group detected"
    
    # Helper function: Check if an atom (except the acid carbon) is attached via a double bond to an oxygen.
    # That would mark it as a carbonyl if not part of the acid group.
    def has_extra_carbonyl(atom, exclude_bond_idxs=set()):
        for bond in atom.GetBonds():
            # Skip bonds already accounted for (e.g. acid carbonyl bond).
            if bond.GetIdx() in exclude_bond_idxs:
                continue
            # Check if bond is double and the other atom is oxygen.
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    # Helper function: Follow a linear chain starting from start_atom (neighbor to acid carbon).
    # We only follow C–C bonds, require no branching, and disallow atoms that are in rings.
    def traverse_linear_chain(start_atom, parent_idx):
        chain_length = 1  # counts start_atom itself
        current_atom = start_atom
        prev_atom_idx = parent_idx

        while True:
            # At current atom, check for extraneous carbonyl functionality.
            # If the atom itself is double-bonded to an oxygen (other than via the acid pattern) then disqualify.
            if has_extra_carbonyl(current_atom):
                return None  # invalid chain candidate

            # Get neighboring carbons (exclude the atom we just came from)
            nbrs = [nbr for nbr in current_atom.GetNeighbors() 
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom_idx]
            # Also, if any atom in the candidate chain is in a ring, disqualify this candidate.
            if current_atom.IsInRing():
                return None

            if len(nbrs) == 0:
                # End of chain reached
                break
            elif len(nbrs) > 1:
                # Branching detected: not unbranched.
                return None
            else:
                # Exactly one neighbor; continue the chain.
                next_atom = nbrs[0]
                prev_atom_idx = current_atom.GetIdx()
                current_atom = next_atom
                chain_length += 1
        return chain_length

    # For each free acid candidate, try to find a valid candidate chain.
    # In our acid match from the SMARTS, match[0] is the carbonyl carbon (the one in [CX3](=O)[OX2H]).
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

        # Identify candidate chain starting atoms: these are carbon neighbors of the acid carbon.
        candidate_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not candidate_neighbors:
            continue  # no attached carbon chain; try next acid candidate

        for start_atom in candidate_neighbors:
            # We only want to traverse candidates that are not in a ring.
            if start_atom.IsInRing():
                continue
            chain_len = traverse_linear_chain(start_atom, acid_carbon_idx)
            if chain_len is None:
                continue  # this candidate chain is disqualified (due to branching, ring involvement or extra carbonyl)
            # total chain length includes the acid carbon.
            total_chain_length = chain_len + 1
            if 4 <= total_chain_length <= 28:
                reason = (f"Found free acid group (–C(=O)[OH]) with an unbranched contiguous aliphatic chain of "
                          f"{total_chain_length} carbons (including the acid carbon).")
                return True, reason
            # Otherwise, this candidate chain does not meet the length criteria—try next candidate.

    return False, "No candidate free acid group is attached to an acceptable unbranched aliphatic chain (4-28 carbons)"

# Example test cases (executed when running this script)
if __name__ == '__main__':
    test_set = [
        # True positives (expected fatty acids)
        ("OC(=O)CCCCCCCCCCC=C", "12-tridecenoic acid"),
        ("C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O", "aspirin-triggered resolvin D2"),
        ("O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O", "Avenoleic acid"),
        ("OC(=O)C=CC=CC=CCCCCCCCCCCCCCCCCC", "Tetracosatrienoic acid"),
        ("C(CCCCCCC[C@H]([C@H](CCCCCCCC)O)O)(=O)O", "(9R,10S)-dihydroxyoctadecanoic acid"),
        ("CCCCCC\\C=C/CCCCCC(O)=O", "cis-tetradec-7-enoic acid"),
        ("OC(=O)CCCCCCCCCCC/C=C/CC", "13-hexadecenoic acid"),
        ("C(O)(=O)CCCCCCCCC(CCCCCCCCC)=O", "10-oxo-nonadecanoic acid"),
        ("CCCCCCCCCC(O)CCC(O)=O", "4-hydroxylauric acid"),
        ("OC(=O)CCC#C/C=C\\C=C\\CCCCCCCCCC#CCCC", "6Z,8E-tricosdien-4,19-diynoic acid"),
        ("CCC=CCC=CCC=CCC=CCCC(=O)O", "4,7,10,13-hexadecatetraenoic acid"),
        ("OC(=O)CCCCCCC/C=C\\CCCCCCCCCC/C=C\\CCCCCC", "28:2(9Z,21Z)"),
        ("C(/C=C/C=C\\C=C\\[C@@H](OO)CCCCC)=C\\[C@H]([C@H](CCCC(O)=O)O)O", "(5S,6R)-dihydroxy-(15S)-hydroperoxy-(7E,9E,11Z,13E)-icosatetraenoic acid"),
        ("C(C(O)=O)C/C=C\\C/C=C\\C/C=C\\C\\C=C/C=C/C(C/C=C\\CC)=O", "(4Z,7Z,10Z,13Z,15E,19Z)-17-oxodocosahexaenoic acid"),
        ("OC(CCCCCC(O)=O)C#CCCCCCCCC", "7-hydroxy-10-heptadecen-8-ynoic acid"),
        ("CCCCCC[C@@H](O)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O", "15(R)-HETE"),
        ("[H]\\C(C)=C\\C(O)=O", "isocrotonic acid"),
        # False positives (should be rejected)
        ("OC(=O)\\C=C\\C(O)=C\\C=O", "(2E,4Z)-4-hydroxy-6-oxohexa-2,4-dienoic acid"),
        ("CCCCCC/C=C\\C/C=C\\CCCCCCCCCCCCC(O)=O", "Prenateic acid"),
        ("OC(CCCCCCCCCC/C=C\\CCC)=O", "Criegeenic acid"),
        ("O[C@@H](CCCC(O)=O)\\C=C/C=C/C=C/[C@H](O)C\\C=C/CCC(O)=O", "18-hydroxy-18-oxo-dinorleukotriene B4"),
        ("OC[C@@H](C(C(O)=O)=O)O", "(S)-3,4-dihydroxy-2-oxobutanoic acid"),
    ]
    
    for smi, name in test_set:
        result, explanation = is_fatty_acid(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nExplanation: {explanation}\n")