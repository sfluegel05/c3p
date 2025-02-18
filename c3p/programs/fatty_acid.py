"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: fatty acid
Definition: “Any aliphatic monocarboxylic acid derived from or contained in esterified form 
in an animal or vegetable fat, oil or wax. Natural fatty acids commonly have a chain of 4 to 28 carbons 
(usually unbranched and even-numbered), which may be saturated or unsaturated. By extension, the term 
is sometimes used to embrace all acyclic aliphatic carboxylic acids.”
This program uses a few heuristic rules:
  1. The molecule must be valid.
  2. It must be “aliphatic” (i.e. contain no rings).
  3. It must contain a carboxyl group. We try to detect:
       - a free acid [C(=O)[OH]]
       - an anionic acid [C(=O)[O-]]
       - or an “ester-like” carbonyl [C(=O)O] that can be part of a fatty acyl chain.
  4. Starting at the carbonyl carbon of the candidate acid group, we “walk” along a contiguous chain 
     of carbon atoms (ignoring the carboxyl oxygen) and require that the chain length is at least 4.
If a candidate is found that meets these criteria, we return True and a brief reason.
Otherwise, the function returns False with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Classifies a molecule as a fatty acid if it meets some basic structural criteria.
    
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      bool: True if the molecule is classified as a fatty acid, otherwise False.
      str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # For our purposes fatty acids are acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and is not acyclic aliphatic"

    # Define SMARTS patterns for candidate acid groups.
    # free acid group: -C(=O)OH
    free_acid = Chem.MolFromSmarts("C(=O)[OH]")
    # anionic acid: -C(=O)[O-]
    anionic_acid = Chem.MolFromSmarts("C(=O)[O-]")
    # ester carbonyl: here we assume that if an ester is present, the acyl part may be the fatty acid.
    ester_carbonyl = Chem.MolFromSmarts("C(=O)O")
    
    # Gather candidate matches. For each match, we note the carbonyl carbon index (first atom in the pattern).
    candidates = []
    for pattern, kind in [(free_acid, "free acid"), (anionic_acid, "anionic acid"), (ester_carbonyl, "ester carbonyl")]:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # In our SMARTS the first atom is the carbonyl carbon.
            candidates.append( (match[0], kind) )

    if not candidates:
        return False, "No carboxyl or ester carbonyl group detected"

    # We now define a helper to compute the longest contiguous chain of carbon atoms
    # starting from a given carbon atom. (The acid carbon is considered part of the chain.)
    def longest_carbon_chain(atom, visited):
        max_length = 1  # count the starting atom itself
        for nbr in atom.GetNeighbors():
            # Only follow carbon atoms
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                new_visited = visited | {nbr.GetIdx()}
                length = 1 + longest_carbon_chain(nbr, new_visited)
                if length > max_length:
                    max_length = length
        return max_length

    # Now, for each candidate acid (or ester) detected, look for a long enough aliphatic chain.
    # Most natural fatty acids have at least 4 carbons (including the carbonyl carbon).
    for carbon_idx, candidate_type in candidates:
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        chain_length = longest_carbon_chain(carbon_atom, {carbon_idx})
        if chain_length < 4:
            # This candidate does not provide a long enough chain.
            continue
        # Optionally, you might want to require a maximum chain length if desired.
        # For our extended definition we accept any chain length >=4.
        reason = (f"Found candidate {candidate_type} with a contiguous aliphatic chain of {chain_length} carbons "
                  f"starting at the carbonyl carbon.")
        return True, reason

    # If we iterated over all candidates without success:
    return False, "No candidate acid group is attached to a long enough aliphatic chain (min 4 carbons)"

# For testing purposes, you might uncomment the code below:
# if __name__ == '__main__':
#     test_smiles = [
#         "OC(=O)CCCCCCCCCCC=C",  # 12-tridecenoic acid, should be fatty acid.
#         "Cl.OC(=O)CCCN1CCC(=CC1)C2=CC=CC=C2",  # contains aromatic ring -> not fatty acid
#         "C1CCCCC1C(=O)O"  # cyclohexanecarboxylic acid, not aliphatic acyclic.
#     ]
#     for smi in test_smiles:
#         result, explanation = is_fatty_acid(smi)
#         print(f"SMILES: {smi}\nResult: {result}\nExplanation: {explanation}\n")