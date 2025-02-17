"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide – any lipid carrying one or more hydroperoxy (-OOH) substituents.
Heuristics (improved):
  1. Look for a hydroperoxy group using a SMARTS pattern that requires an sp3 carbon attached to an O–O–H,
     i.e. [CX4]-[OX2]-[OX2H]. Also check that each atom in the group AND its immediate neighbors (within one bond)
     are formally neutral. This helps to avoid matching groups where part of a carboxylate is deprotonated.
  2. From the carbon which carries the −OOH group, perform a DFS to find a contiguous aliphatic (non‐ring, non‐aromatic)
     chain. The required minimum chain length is lowered to 6 carbons.
  3. Finally, require that the molecular weight is above 250 Da.
If all conditions are met, the molecule is classified as a lipid hydroperoxide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule represents a lipid hydroperoxide.
    A lipid hydroperoxide here is defined as a molecule that carries at least one neutral hydroperoxy (-OOH)
    substituent attached to a sufficiently long (>=6 carbons) contiguous, aliphatic (non-aromatic, non-ring) segment,
    and has a molecular weight above a minimal threshold.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lipid hydroperoxide; False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Look for a hydroperoxy group using SMARTS.
    # We require an sp3 carbon attached to an O–O–H group.
    hydroperoxy_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[OX2H]")
    hydro_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    valid_hydro_match = None
    
    # Helper function: check that an atom and its immediate neighbors are all neutral.
    def local_environment_is_neutral(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetFormalCharge() != 0:
            return False
        for nb in atom.GetNeighbors():
            if nb.GetFormalCharge() != 0:
                return False
        return True

    # Try to find a matching hydroperoxy group meeting our local neutrality criteria.
    for match in hydro_matches:
        # match is a tuple: (hydro_carbon, middle oxygen, terminal oxygen with H)
        if (local_environment_is_neutral(match[0]) and 
            local_environment_is_neutral(match[1]) and 
            local_environment_is_neutral(match[2])):
            valid_hydro_match = match
            break
    if valid_hydro_match is None:
        return False, "No valid neutral hydroperoxy (-OOH) group found"
    
    # 3. Check that there is a sufficiently long contiguous aliphatic chain.
    # Starting from the carbon atom which carries the hydroperoxy group (first atom of the match),
    # we perform a DFS. We allow only carbon atoms that are non-aromatic and not in any ring.
    # The required length is lowered to 6 (so that we do not miss some lipids such as prostaglandins).
    def longest_chain_depth(atom_idx, visited):
        """Recursively compute the maximum chain length from a given carbon atom.
           Only traverse to neighbors that are carbons, non-aromatic and not in a ring.
           Returns the chain length counting the starting atom.
        """
        visited.add(atom_idx)
        current_atom = mol.GetAtomWithIdx(atom_idx)
        max_depth = 1  # count the current atom
        # Explore each neighboring carbon.
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

    hydro_carbon_idx = valid_hydro_match[0]
    chain_length = longest_chain_depth(hydro_carbon_idx, set())
    if chain_length < 6:
        return False, f"Longest contiguous aliphatic chain (length {chain_length}) attached to the hydroperoxy group is too short"
    
    # 4. Check overall molecular weight; many lipids have MW > 250 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical lipid"
    
    return True, "Molecule is classified as a lipid hydroperoxide: contains a neutral hydroperoxy (-OOH) group attached to a long aliphatic chain"

# Example usage:
if __name__ == "__main__":
    # Test with one of the example molecules.
    test_smiles = "O(O)[C@@H](CCCC(O)=O)/C=C/C=C\\C/C=C\\C/C=C\\CC"  # 5S-HpEPE
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)