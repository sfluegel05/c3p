"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: fatty acid
Definition: “Any acyclic aliphatic carboxylic acid. 
Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched and even‐numbered), 
which may be saturated or unsaturated. By extension, the term is sometimes used to embrace all acyclic aliphatic carboxylic acids.”
This classifier rejects molecules that are cyclic or that contain any nitrogen (to avoid peptides/amino acid derivatives) 
and then requires that a free (undissociated) carboxyl group –C(=O)OH is attached to a contiguous, unbranched carbon chain of 4 to 28 carbons.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if the given SMILES is a fatty acid.
    For this classifier, a fatty acid must be:
      - A valid molecule,
      - Acyclic (no rings) and free of nitrogen atoms,
      - Contain a free (undissociated) carboxylic acid group (–C(=O)OH),
      - Have an unbranched, contiguous aliphatic chain (starting at the carbonyl carbon) whose total length is between 4 and 28 carbons.
    
    Args:
      smiles (str): SMILES string for the molecule.
    
    Returns:
      bool: True if the molecule is classified as a fatty acid, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Rule 1: For our fatty acid definition, molecule must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and is not acyclic aliphatic"
    
    # Rule 2: Do not allow nitrogen atoms (indicative of peptides or amino acids)
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Molecule contains nitrogen atoms indicating peptide or amino acid derivatives"
    
    # Rule 3: Look for free carboxylic acid group (–C(=O)OH).
    # SMARTS: [CX3](=O)[OX2H] matches a free (undissociated) carboxyl.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No free carboxyl (–C(=O)OH) group detected"
    
    # Helper function: Recursively compute the length of a linear, unbranched carbon chain.
    # The chain should only follow carbon neighbors that are not in the current path.
    def linear_chain_length(atom, parent_idx, visited):
        # Count the current atom as 1.
        length = 1
        # List carbon neighbors excluding the one we just came from.
        nbrs = [nbr for nbr in atom.GetNeighbors() 
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != parent_idx]
        if len(nbrs) == 0:
            return length
        # If more than one carbon neighbor appears, then the chain is branched.
        if len(nbrs) > 1:
            return None  # Indicates branching
        # Continue the chain if not visited.
        next_atom = nbrs[0]
        if next_atom.GetIdx() in visited:
            return length
        visited.add(next_atom.GetIdx())
        cont = linear_chain_length(next_atom, atom.GetIdx(), visited)
        if cont is None:
            return None  # Propagate failure due to branching further down
        return length + cont

    # For each free acid candidate, attempt to traverse a linear chain starting at the carbonyl carbon.
    # In our acid match, match[0] is the carbonyl carbon.
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        
        # Get carbon neighbors of the acid carbon (exclude oxygens)
        candidate_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not candidate_neighbors:
            continue  # no alkyl chain attached to acid carbon
        
        # We then consider each neighboring carbon as the starting point of the contiguous chain.
        for start_atom in candidate_neighbors:
            visited = {acid_carbon_idx, start_atom.GetIdx()}  # start atom visited along with acid carbon
            chain_len = linear_chain_length(start_atom, acid_carbon_idx, visited)
            if chain_len is None:
                # This candidate chain is branched.
                continue
            # We include the acid (carbonyl) carbon in the count.
            total_chain = chain_len + 1
            # Check if the chain length is in the acceptable range.
            if 4 <= total_chain <= 28:
                reason = (f"Found free acid group (–C(=O)OH) with an unbranched contiguous aliphatic chain of "
                          f"{total_chain} carbons (including the acid carbon).")
                return True, reason
            # If the chain is too short or too long, report that candidate but continue trying others.
    return False, "No candidate free acid group is attached to an acceptable unbranched aliphatic chain (4-28 carbons)"

# Example tests (can be run when executing the script)
if __name__ == '__main__':
    test_set = [
        # True positives (expected fatty acids)
        ("OC(=O)CCCCCCCCCCC=C", "12-tridecenoic acid"),
        ("C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O", "aspirin-triggered resolvin D2"),
        ("O=C(O)[C@H]([C@@H](C(=O)O)C)CCCCCCCC", "Sphaeric acid"),
        ("O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O", "Avenoleic acid"),
        ("OC(=O)C=CC=CC=CCCCCCCCCCCCCCCCCC", "Tetracosatrienoic acid"),
        ("CC(C)C[C@@H](O)C(O)=O", "(R)-2-hydroxy-4-methylpentanoic acid"),
        ("CCCCCC\\C=C/CCCCCC(O)=O", "cis-tetradec-7-enoic acid"),
        ("OC(=O)CCCCCCCCCCC/C=C/CC", "13-hexadecenoic acid"),
        ("C(O)(=O)CCCCCCCCC(CCCCCCCCC)=O", "10-oxo-nonadecanoic acid"),
        ("CCCCCCCCCC(O)CCC(O)=O", "4-hydroxylauric acid"),
        ("OC(=O)CCC(CCCC)CC", "4-Ethyloctanoic acid"),
        ("O=C(CCCCCCCC(O)=O)/C=C/C=C\\C/C=C\\CC", "9-OxoOTrE"),
        ("OC(=O)CCC#C/C=C\\C=C\\CCCCCCCCCC#CCCC", "6Z,8E-tricosdien-4,19-diynoic acid"),
        ("CCC=CCC=CCC=CCC=CCCC(=O)O", "4,7,10,13-hexadecatetraenoic acid"),
        ("OC(=O)CCCCCCC/C=C\\CCCCCCCCCC/C=C\\CCCCCC", "28:2(9Z,21Z)"),
        ("C(/C=C/C=C\\C=C\\[C@@H](OO)CCCCC)=C\\[C@H]([C@H](CCCC(O)=O)O)O", "(5S,6R)-dihydroxy-(15S)-hydroperoxy-(7E,9E,11Z,13E)-icosatetraenoic acid"),
        ("C(C(O)=O)C/C=C\\C/C=C\\C/C=C\\C\\C=C/C=C/C(C/C=C\\CC)=O", "(4Z,7Z,10Z,13Z,15E,19Z)-17-oxodocosahexaenoic acid"),
        ("OC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "Hentriacontanoic acid"),
        ("OC(=O)[C@H](C[C@H](CCCCCCCCCCCCCCCCCC)C)C", "Mycosanoic acid (C24)"),
        ("OC(CCCCCC(O)=O)C#CCCCCCCCC", "7-hydroxy-10-heptadecen-8-ynoic acid"),
        ("CCCCCC[C@@H](O)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O", "15(R)-HETE"),
        ("[H]\\C(C)=C\\C(O)=O", "isocrotonic acid"),
        # False positives (should be rejected)
        ("OC1C(C(C(=O)C1)C/C=C\\C/C=C\\CC)/C=C/C(O)C/C=C\\CCC(O)=O", "7-hydroxy-D4-neuroprostane"),
        ("CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O", "juvenile hormone I acid"),
        ("O=C(O)/C=C/C=C/C=C/C=C/[C@@]1([C@@H](C(=C[C@H]2[C@H]1CC[C@@H]([C@@H]2O)C)CO)/C(=C/C)/C)C", "Janthinopolyenemycin B"),
        ("O1C(C1)(C(OC(=O)C(O)C(CC)C)/C=C/C=C/C=C/C(O)=O)C", "AF Toxin II"),
        ("Cl.OC(=O)CCCN1CCC(=CC1)C2=CC=CC=C2", "4-[4-phenyl-3,6-dihydro-1(2H)-pyridinyl]butanoic acid hydrochloride"),
        ("OC(CCC(C)C)CC(O)=O", "3-hydroxy-isoheptanoic acid"),
        ("OC(=O)\\C=C\\C(O)=C\\C=O", "(2E,4Z)-4-hydroxy-6-oxohexa-2,4-dienoic acid"),
        ("CCCCCC/C=C\\C/C=C\\CCCCCCCCCCCCC(O)=O", "Prenateic acid"),
    ]
    
    for smi, name in test_set:
        result, explanation = is_fatty_acid(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nExplanation: {explanation}\n")