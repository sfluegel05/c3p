"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: fatty acid
Definition: “Any acyclic aliphatic carboxylic acid. Natural fatty acids commonly have a chain of 4 to 28 carbons,
which may be saturated or unsaturated.”
This improved version uses refined rules:
  1. The molecule must be valid and acyclic.
  2. It must contain a free (undissociated) carboxyl group (–C(=O)OH) rather than an ester or anionic acid.
  3. Starting from the carbonyl carbon of that group, we walk along a contiguous chain of carbon atoms (ignoring
     the carboxyl oxygen). The chain must be at least 4 carbons long.
If these criteria are met, it is classified as a fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    For this classifier, a fatty acid must be an acyclic molecule that contains a free carboxylic acid group
    (–C(=O)OH) and a contiguous aliphatic chain (starting at the carbonyl carbon) of at least 4 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a fatty acid, False otherwise.
        str: A brief explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Rule 1: Fatty acids are acyclic aliphatic molecules.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and is not acyclic aliphatic"

    # Rule 2: Search for a free (undissociated) carboxylic acid group.
    # We define a SMARTS that matches a carbonyl carbon [CX3] double-bonded to an oxygen and single bonded to
    # an oxygen that has at least one hydrogen (i.e. –OH). In many cases the hydrogen is implicit.
    free_acid_smarts = "[CX3](=O)[OX2H]" 
    free_acid = Chem.MolFromSmarts(free_acid_smarts)
    acid_matches = mol.GetSubstructMatches(free_acid)
    if not acid_matches:
        return False, "No free carboxyl (–C(=O)OH) group detected"

    # Helper: From a given starting carbon atom, compute the length of the longest contiguous chain of carbons.
    def longest_carbon_chain(atom, visited):
        max_length = 1  # count the starting atom
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                # Continue along any carbon, regardless of bond order (saturated or unsaturated)
                new_visited = visited | {nbr.GetIdx()}
                chain_length = 1 + longest_carbon_chain(nbr, new_visited)
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # Rule 3: For each acid candidate, compute the contiguous chain length starting at the carbonyl carbon.
    # Our SMARTS pattern is arranged so that the first atom in the match is the carbonyl carbon.
    for match in acid_matches:
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        chain_length = longest_carbon_chain(carbon_atom, {carbon_idx})
        if chain_length >= 4:
            reason = (f"Found free acid group (–C(=O)OH) with a contiguous aliphatic chain of "
                      f"{chain_length} carbons starting at the carbonyl carbon.")
            return True, reason

    return False, "No candidate free acid group is attached to a long enough aliphatic chain (min 4 carbons)"

# Example tests:
if __name__ == '__main__':
    test_set = [
        # True positives
        ("OC(=O)CCCCCCCCCCC=C", "12-tridecenoic acid"),
        ("C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O", "aspirin-triggered resolvin D2"),
        ("O=C(O)[C@H]([C@@H](C(=O)O)C)CCCCCCCC", "Sphaeric acid"),
        ("O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O", "Avenoleic acid"),
        ("OC(=O)C=CC=CC=CCCCCCCCCCCCCCCCCC", "Tetracosatrienoic acid"),
        ("CC(C)C[C@@H](O)C(O)=O", "(R)-2-hydroxy-4-methylpentanoic acid"),
        ("CCCCCC\\C=C/CCCCCC(O)=O", "cis-tetradec-7-enoic acid"),
        ("OC(=O)CCCCCCCCCCC/C=C/CC", "13-hexadecenoic acid"),
        ("C(O)(=O)CCCCCCCCC(CCCCCCCCC)=O", "10-oxo-nonadecanoic acid"),
        ("OC(=O)CCCN(N=O)CCC", "4-[nitroso(propyl)amino]butanoic acid"),
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
        # False due to rings
        ("OC1C(C(C(=O)C1)C/C=C\\C/C=C\\CC)/C=C/C(O)C/C=C\\CCC(O)=O", "7-hydroxy-D4-neuroprostane"),
        ("CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O", "juvenile hormone I acid"),
    ]
    
    for smi, name in test_set:
        result, explanation = is_fatty_acid(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nExplanation: {explanation}\n")