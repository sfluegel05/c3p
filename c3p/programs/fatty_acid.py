"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: fatty acid
Definition: “Any acyclic aliphatic carboxylic acid. 
Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched), which may be saturated or unsaturated.”
The classifier now rejects molecules that are cyclic or that contain any nitrogen atoms (e.g. peptides).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    For this classifier, a fatty acid must be:
      - A valid molecule,
      - Acyclic (no rings) and free of nitrogen atoms (to avoid peptides or amino-acid derivatives),
      - Have a free (undissociated) carboxylic acid group (–C(=O)OH),
      - Connected to a contiguous aliphatic chain (starting at the carbonyl carbon) of at least 4 carbons.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as a fatty acid, False otherwise.
        str: An explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Rule 1: Fatty acids are acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and is not acyclic aliphatic"
        
    # Extra check: Do not allow nitrogen atoms (which usually indicate peptides or amino acid derivatives).
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Molecule contains nitrogen atoms indicating peptide or amino acid derivatives"
    
    # Rule 2: The molecule must contain a free carboxylic acid group.
    # SMARTS here covers a free (undissociated) carboxyl: carbonyl carbon (CX3) bonded to an –OH (OX2H).
    free_acid_smarts = "[CX3](=O)[OX2H]"
    free_acid = Chem.MolFromSmarts(free_acid_smarts)
    acid_matches = mol.GetSubstructMatches(free_acid)
    if not acid_matches:
        return False, "No free carboxyl (–C(=O)OH) group detected"
    
    # Helper function: From a given starting carbon atom, compute the length of the longest contiguous carbon chain.
    def longest_carbon_chain(atom, visited):
        max_length = 1  # counting the starting atom
        for nbr in atom.GetNeighbors():
            # Only follow carbon atoms that haven’t been visited
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                new_visited = visited | {nbr.GetIdx()}
                chain_length = 1 + longest_carbon_chain(nbr, new_visited)
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # Rule 3: For each acid candidate, compute the contiguous chain length from the carbonyl carbon.
    # With our SMARTS the first atom in the match is the carbonyl carbon.
    for match in acid_matches:
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        chain_length = longest_carbon_chain(carbon_atom, {carbon_idx})
        if chain_length >= 4:
            reason = (f"Found free acid group (–C(=O)OH) with a contiguous aliphatic chain of "
                      f"{chain_length} carbons starting at the carbonyl carbon.")
            return True, reason

    return False, "No candidate free acid group is attached to a long enough (>=4 carbons) aliphatic chain"

# Example tests:
if __name__ == '__main__':
    test_set = [
        # True positives (fatty acids)
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
        # False positives (should be rejected)
        ("OC1C(C(C(=O)C1)C/C=C\\C/C=C\\CC)/C=C/C(O)C/C=C\\CCC(O)=O", "7-hydroxy-D4-neuroprostane"),
        ("CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O", "juvenile hormone I acid"),
        ("O=C(O)/C=C/C=C/C=C/C=C/[C@@]1([C@@H](C(=C[C@H]2[C@H]1CC[C@@H]([C@@H]2O)C)CO)/C(=C/C)/C)C", "Janthinopolyenemycin B"),
        ("O1C(C1)(C(OC(=O)C(O)C(CC)C)/C=C/C=C/C=C/C(O)=O)C", "AF Toxin II"),
        ("Cl.OC(=O)CCCN1CCC(=CC1)C2=CC=CC=C2", "4-[4-phenyl-3,6-dihydro-1(2H)-pyridinyl]butanoic acid hydrochloride"),
        ("CCCCCC/C=C\\C/C=C\\CCCCCCCCCCCCC(O)=O", "Prenateic acid"),
    ]
    
    for smi, name in test_set:
        result, explanation = is_fatty_acid(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nExplanation: {explanation}\n")