"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide – any lipid carrying one or more hydroperoxy (-OOH) substituents.
Improved heuristics:
  1. Parse the SMILES and reject molecules that are too small (<250 Da) or not neutral.
  2. Look for at least one neutral hydroperoxy group using a SMARTS pattern on an sp3 carbon.
  3. Look for a carboxyl/ester functional group as a marker of a fatty acid moiety.
  4. For each such group, traverse from the carbonyl carbon over non‐aromatic, non‐ring (aliphatic) carbons
     to find the longest contiguous acyl chain. A minimum length (>=8 carbons) is required.
  5. If all conditions are met, the molecule is classified as a lipid hydroperoxide.
  
In our evaluation the error in the last attempt was likely because the hydroperoxy group check was too local;
thus we now “anchor” our search by requiring the presence of a fatty acid moiety.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule represents a lipid hydroperoxide.
    A lipid hydroperoxide is defined here as a molecule that:
      - is uncharged and has a reasonable molecular weight (>250 Da);
      - contains at least one neutral hydroperoxy (-OOH) group 
        attached to an sp3 (non‐aromatic) carbon (via pattern matching);
      - also carries a fatty acid (or ester) moiety as evidenced by 
        the presence of a carboxyl/ester group whose acyl chain (traversed over non‐aromatic, non‐ring carbons)
        is at least 8 carbons long.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lipid hydroperoxide; False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Parse SMILES and basic checks.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Ensure the molecule overall is uncharged.
    if Chem.GetFormalCharge(mol) != 0:
        return False, f"Molecule has a net formal charge of {Chem.GetFormalCharge(mol)}, expected 0"
    # Check for a minimum molecular weight – many lipids are large (>250 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical lipid"
    
    # 2. Find a hydroperoxy group.
    # We require a neutral hydroperoxy group: using a SMARTS for an sp3 carbon (non‐aromatic) attached to an –OOH.
    hydroperoxy_pattern = Chem.MolFromSmarts("[#6;!a]-[OX2]-[OX2H]")
    hydro_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    
    def local_environment_is_neutral(atom_idx):
        """Check that the atom and its immediate neighbors are all formally neutral."""
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetFormalCharge() != 0:
            return False
        for nb in atom.GetNeighbors():
            if nb.GetFormalCharge() != 0:
                return False
        return True
    
    valid_hydro_match = None
    for match in hydro_matches:
        # match: (hydrocarbon, middle oxygen, terminal oxygen with H)
        if (local_environment_is_neutral(match[0]) and 
            local_environment_is_neutral(match[1]) and 
            local_environment_is_neutral(match[2])):
            valid_hydro_match = match
            break
    if valid_hydro_match is None:
        return False, "No valid neutral hydroperoxy (-OOH) group found"
    
    # 3. Look for a carboxyl/ester group as a marker for a fatty acid moiety.
    # We use SMARTS that covers carboxylic acids and esters: 
    #   [CX3](=O)[OX1H] covers acids while [CX3](=O)[OX2H0] covers esters.
    carboxyl_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1H,OX2H0]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_smarts)
    if not carboxyl_matches:
        return False, "No carboxyl or ester functional group found; not a typical lipid structure"
    
    # 4. Check that at least one of the carboxyl/ester groups is connected
    #    to a sufficiently long contiguous aliphatic (non-aromatic, non-ring) chain.
    # We define a DFS that only traverses carbon atoms meeting these conditions.
    def longest_chain_depth(atom_idx, visited):
        """Recursively computes the maximum contiguous chain length from a given atom,
           traversing only carbons that are non-aromatic and not in rings.
           Returns the count including the starting atom.
        """
        visited.add(atom_idx)
        current_atom = mol.GetAtomWithIdx(atom_idx)
        max_depth = 1  # count the current atom
        for neighbor in current_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if (neighbor.GetAtomicNum() == 6 and
                not neighbor.GetIsAromatic() and
                not neighbor.IsInRing() and
                n_idx not in visited):
                depth = 1 + longest_chain_depth(n_idx, visited.copy())
                if depth > max_depth:
                    max_depth = depth
        return max_depth
    
    chain_lengths = []
    # For each carboxyl/ester match, the first atom in the match is the carbonyl carbon.
    for match in carboxyl_matches:
        carbonyl_idx = match[0]
        chain_length = longest_chain_depth(carbonyl_idx, set())
        chain_lengths.append(chain_length)
    
    if not chain_lengths:
        return False, "No contiguous aliphatic chain found from any carboxyl group"
    
    max_chain = max(chain_lengths)
    if max_chain < 8:
        return False, f"Longest contiguous acyl chain is too short (length {max_chain}); expected >= 8 carbons"
    
    return True, ("Molecule is classified as a lipid hydroperoxide: contains a neutral hydroperoxy (-OOH) group "
                  "and a fatty acid (or ester) moiety with a long aliphatic chain.")

# Example usage:
if __name__ == "__main__":
    # Try one example: 9-HpOTrE (should be classified as lipid hydroperoxide).
    test_smiles = "O(O)C(CCCCCCCC(O)=O)/C=C/C=C/C/C=C\\CC"
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)