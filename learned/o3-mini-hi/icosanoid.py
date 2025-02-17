"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation of the three C20 essential fatty acids (EFAs) – icosapentaenoic acid (EPA), arachidonic acid (AA) and dihomo-gamma-linolenic acid (DGLA).
Heuristic criteria applied:
  • The SMILES must be valid.
  • The overall molecular weight should lie between 250 and 900 Da.
  • The total carbon count must be at least 15 (but we do not enforce a strict upper limit, since conjugation may add extra carbons).
  • The molecule must contain at least one oxygenated carbonyl fragment (as in an acid or ester group).
  • It must contain at least one non-aromatic C=C double bond (to help avoid counting carbonyl C=O bonds or aromatic systems).
  • Its longest continuous carbon chain (as a proxy for the oxidized fatty acid backbone) should be roughly between 12 and 25 carbons.
Note: This is a heuristic that may miss nuances of icosanoid chemistry.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the longest continuous chain containing only carbon atoms.
    Uses a depth-first search (DFS) on the subgraph of carbon atoms.
    """
    # Build graph: for every carbon, list its carbon neighbors.
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            idx = atom.GetIdx()
            carbon_graph[idx] = []
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    carbon_graph[idx].append(nbr.GetIdx())
    
    # DFS to find the longest simple path in the carbon graph.
    max_length = 0
    def dfs(node, visited):
        current_max = 1
        for neighbor in carbon_graph.get(node, []):
            if neighbor not in visited:
                new_length = 1 + dfs(neighbor, visited | {neighbor})
                if new_length > current_max:
                    current_max = new_length
        return current_max

    for node in carbon_graph:
        chain_length = dfs(node, {node})
        if chain_length > max_length:
            max_length = chain_length
    return max_length
    
    
def is_icosanoid(smiles: str):
    """
    Determines if a molecule could be an icosanoid based on its SMILES string.
    
    Heuristic filters:
      - Valid molecule.
      - Molecular weight between 250 and 900 Da.
      - Total carbon atom count of at least 15.
      - Contains a carbonyl-oxygen fragment (as in carboxylic acid/ester, using SMARTS "C(=O)O").
      - Contains at least one non-aromatic C=C (alkene) double bond.
      - Has a longest continuous carbon chain (computed via DFS) between 12 and 25 atoms long.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule passes our heuristic for being an icosanoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count total carbon atoms.
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_c < 15:
        return False, f"Total carbon count {total_c} is too low for an icosanoid."
    
    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is out of expected 250-900 Da range."
    
    # Check for oxygenated carbonyl fragment using SMARTS "C(=O)O".
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl-oxygen fragment (acid/ester group) found."
    
    # Count non-aromatic C=C double bonds.
    alkene_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            # Check that both atoms are carbon and that the bond is not aromatic.
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and not bond.GetIsAromatic():
                alkene_count += 1
    if alkene_count < 1:
        return False, "Too few non-aromatic C=C double bonds; icosanoids are typically polyunsaturated."
    
    # Compute the longest continuous carbon chain.
    chain_length = longest_carbon_chain(mol)
    if chain_length < 12 or chain_length > 25:
        return False, f"Longest continuous carbon chain length ({chain_length}) is not in expected range (12-25) for an icosanoid."
    
    # Check for free acid group (free carboxyl) which many icosanoids have.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        acid_note = " (No free carboxylic acid detected; molecule may be fully esterified/conjugated.)"
    else:
        acid_note = ""
    
    return True, (f"Meets heuristic criteria for an icosanoid: "
                  f"(Total C={total_c}, MW={mol_wt:.1f} Da, non-aromatic C=C count={alkene_count}, "
                  f"Longest C-chain={chain_length}){acid_note}")


# Example usage:
if __name__ == "__main__":
    # Some examples from the provided list.
    smiles_examples = [
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CCCC[C@H](OO)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O",      # 15(S)-HPETE
        "O[C@H]1[C@@H]([C@H]([C@@H](O)C1)C/C=C\\CCCC(OC(C)C)=O)CC[C@@H](O)CCOC2=CC=CC=C2",  # 17-phenoxy trinor PGF2alpha isopropyl ester
        "[H]C(CC)=C([H])CCCCCCCCCCCCCCCC(O)=O"  # 17-icosenoic acid
    ]
    
    for smi in smiles_examples:
        result, reason = is_icosanoid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")