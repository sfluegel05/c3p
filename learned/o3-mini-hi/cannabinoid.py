"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid class

Definition:
  "A diverse group of pharmacologically active secondary metabolite characteristic to Cannabis plant
   as well as produced naturally in the body by humans and animals. Cannabinoids contain oxygen as a 
   part of the heterocyclic ring or in the form of various functional groups."

Heuristic rules in this improved version:
  1. Molecule must be a valid structure, contain oxygen, have no elements (e.g. phosphorus) deemed inappropriate,
     a total carbon count of at least 15, and a molecular weight above ~200 Da.
  2. For phytocannabinoids, the molecule must contain a dihydroxybenzene substructure 
     (SMARTS "c1cc(O)c(O)cc1").
  3. For synthetic cannabinoids, we check for an indole ring using the substructure "c1cc2c(c1)[nH]c(c2)".
  4. For endocannabinoids, if the molecule contains an amide or ester group and a long aliphatic chain
     (here determined by a DFS search for the longest carbon chain, requiring chain length ≥ 16), then it is classified.
  5. As a fallback, if the molecule has a heterocyclic ring (≥5 atoms) that includes oxygen and also has a long carbon chain (≥16),
     then it might be cannabinoid‐related.
  
If none of these rules apply then the molecule is not classified as a cannabinoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Returns the length (number of atoms) of the longest chain 
    of carbon atoms in the molecule. It builds an undirected graph
    of carbon atoms (atomic number==6) and finds the longest simple path.
    """
    # Obtain indices of all carbon atoms.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Build connectivity: only consider bonds between carbons.
    graph = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    
    longest = 0
    # Depth-first search for longest simple path
    def dfs(node, visited, length):
        nonlocal longest
        if length > longest:
            longest = length
        for nbr in graph[node]:
            if nbr not in visited:
                dfs(nbr, visited | {nbr}, length + 1)
                
    for idx in carbon_idxs:
        dfs(idx, {idx}, 1)
    
    return longest

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on selected heuristics.
    
    The function implements these checks:
      - Validity, sufficient oxygen, no phosphorus, a minimum total number of carbons (>=15)
        and molecular weight (>=200 Da).
      - If a dihydroxybenzene substructure (resorcinol motif) is present then the molecule is interpreted
        as a phytocannabinoid.
      - If an indole ring (common in synthetic cannabinoids) is present then it is classified as cannabinoid.
      - Otherwise, if the molecule has an amide or ester group AND the longest carbon chain is at least 16 atoms long,
        then it is treated as an endocannabinoid.
      - As a fallback, if the molecule contains a heterocyclic ring (≥5 atoms) including oxygen and has a long aliphatic
        chain (chain length ≥ 16) then it is flagged as possibly cannabinoid-related.
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cannabinoid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules with phosphorus (atomic number 15)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus which is not typical in cannabinoid structures"
    
    # Ensure molecule contains oxygen.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "Molecule does not contain oxygen, but cannabinoids typically require oxygen."
    
    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, f"Total carbon count ({total_carbons}) is too low to be a typical cannabinoid."
    
    # Check molecular weight (exact weight) is above ~200 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low to be a typical cannabinoid."
    
    # Rule 1: Look for a dihydroxybenzene (resorcinol) motif common in phytocannabinoids.
    dihydroxybenzene = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if dihydroxybenzene and mol.HasSubstructMatch(dihydroxybenzene):
        return True, "Molecule contains a dihydroxybenzene (resorcinol) motif, common in phytocannabinoids."
    
    # Rule 2: Look for an indole ring (common in many synthetic cannabinoids).
    indole = Chem.MolFromSmarts("c1cc2c(c1)[nH]c(c2)")
    if indole and mol.HasSubstructMatch(indole):
        return True, "Molecule contains an indole ring, a feature common in synthetic cannabinoids."
    
    # Rule 3: Check for polar functionality (amide or ester) paired with a long carbon chain.
    amide = Chem.MolFromSmarts("[NX3][CX3](=O)")
    ester = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    has_amide = mol.HasSubstructMatch(amide) if amide else False
    has_ester = mol.HasSubstructMatch(ester) if ester else False
    
    chain_length = longest_carbon_chain(mol)
    
    if (has_amide or has_ester) and chain_length >= 16:
        return True, f"Molecule has a polar head group (amide/ester) and a long aliphatic chain (chain length {chain_length}), features common in endocannabinoids."
    
    # Rule 4: As a fallback, check for any heterocyclic ring (ring of >=5 atoms) containing oxygen
    # along with evidence of a long carbon chain.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) >= 5:
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                if chain_length >= 16:
                    return True, "Molecule contains a heterocyclic ring with oxygen and a long aliphatic chain, a possible cannabinoid feature."
    
    # If none of the rules apply, do not classify as cannabinoid.
    return False, "Molecule does not contain key cannabinoid structural features (dihydroxybenzene or indole core, or polar group with long chain, or oxygenated heterocycle with long chain)."

# Example calls (uncomment to test):
# test_smiles = [
#    "O=C(O[C@@H]([C@@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C",  # Roselipin 3E
#    "O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO",  # MG(0:0/20:2(11Z,14Z)/0:0)
#    "CCCC\\C=C/CCCCCCCC(=O)NCCO",  # palmitoleoyl ethanolamide
#    "CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCOC(CO)CO",  # 2-arachidonyl glyceryl ether
#    "N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF",  # 1-(5-fluoropentyl)-3-(2-iodobenzoyl)indole
# ]
# for smi in test_smiles:
#     result, reason = is_cannabinoid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")