"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: Lipopeptide – a compound consisting of a peptide (i.e. one or more amide bonds)
with an attached lipid (i.e. a long aliphatic chain).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_longest_aliphatic_chain(mol):
    """
    Returns the length of the longest chain (number of atoms) that is:
      – composed solely of nonaromatic carbon atoms (atomic number 6)
      – and where each carbon is sp3-hybridized.
    This DFS walks over the molecular graph (ignoring bonds to heteroatoms).
    """
    longest = 0
    # Cache for DFS results per starting atom (atom idx, coming from neighbor idx)
    memo = {}
    
    def dfs(atom, coming_from):
        key = (atom.GetIdx(), coming_from)
        if key in memo:
            return memo[key]
        
        max_length = 1  # count the current atom
        for neighbor in atom.GetNeighbors():
            # Avoid going back to the atom we came from
            if coming_from is not None and neighbor.GetIdx() == coming_from:
                continue
            # Only follow if this neighbor is a carbon, is sp3 (we check degree and aromaticity),
            # and is not aromatic.
            if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic() and neighbor.GetHybridization().name == "SP3":
                length = 1 + dfs(neighbor, atom.GetIdx())
                if length > max_length:
                    max_length = length
        memo[key] = max_length
        return max_length

    for atom in mol.GetAtoms():
        # Start DFS from each carbon that is nonaromatic, sp3
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and atom.GetHybridization().name == "SP3":
            chain_length = dfs(atom, None)
            if chain_length > longest:
                longest = chain_length
    return longest

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide should contain a peptide component (evidenced by amide bonds)
    and a lipid component (evidenced by a long aliphatic chain attached).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is recognized as a lipopeptide, False otherwise.
        str: A reason detailing the classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, check for peptide (amide) bonds.
    # Typical peptide/amide bond fragment pattern: -N-C(=O)-
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    peptide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(peptide_matches) < 1:
        return False, "No amide (peptide) bonds found; peptide component missing"
    
    # Next, try to detect the lipid component. Rather than a fixed SMARTS,
    # we compute the length of the longest chain of nonaromatic sp3 carbons.
    longest_chain = get_longest_aliphatic_chain(mol)
    if longest_chain < 8:
        return False, f"Longest aliphatic chain has {longest_chain} carbons; lipid component missing"
    
    # If both a peptide signature and a long lipid chain are found, classify as lipopeptide.
    return True, f"Contains amide bonds (peptide component) and a long aliphatic chain of {longest_chain} carbons (lipid component)"

# Example usage (uncomment below lines to test):
# smiles_examples = [
#     "[H][C@@]1(CCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)O1",  # surfactin C
#     "OC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # long fatty acid (should not be classified as lipopeptide)
# ]
# for s in smiles_examples:
#     result, reason = is_lipopeptide(s)
#     print(result, reason)