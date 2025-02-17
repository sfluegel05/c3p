"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: aldopentose, defined as 'A pentose with a (potential) aldehyde group at one end.'
This function checks that the molecule has exactly 5 carbons and 5 oxygens
(consistent with an aldopentose formula, C5H10O5) and then determines whether it is in
the open-chain or cyclized (hemiacetal) form. In the open-chain form, an aldehyde group is
expected at a terminal carbon; in the cyclic form (without lactone functionality) the molecule
is assumed to be in equilibrium with the open‐chain form.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    
    An aldopentose is a 5-carbon sugar (typically C5H10O5) that in its open-chain form displays
    an aldehyde group at one terminal carbon. In its cyclic (hemiacetal) form the aldehyde may not
    be explicit, but the molecule still has 5 carbons in its sugar backbone.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an aldopentose; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count the number of carbon and oxygen atoms in the molecule.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if carbon_count != 5:
        return False, f"Expected 5 carbon atoms for an aldopentose; found {carbon_count}."
    if oxygen_count != 5:
        return False, f"Expected 5 oxygen atoms for an aldopentose; found {oxygen_count}."
    
    # Prepare the aldehyde SMARTS pattern.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    
    # Check ring information.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Define lactone (ester) pattern to exclude cyclic lactones.
    lactone_pattern = Chem.MolFromSmarts("[CX3](=O)[O]")
    
    # For acyclic (open-chain) structures:
    if num_rings == 0:
        # The molecule must have a terminal aldehyde group.
        if not mol.HasSubstructMatch(aldehyde_pattern):
            return False, "Acyclic structure missing an aldehyde group."
        
        # Build a carbon connectivity graph (only for 5 atoms so it's fast).
        carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
        carbon_graph = {idx: set() for idx in carbon_idxs}
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                idx1 = a1.GetIdx()
                idx2 = a2.GetIdx()
                if idx1 in carbon_graph and idx2 in carbon_graph:
                    carbon_graph[idx1].add(idx2)
                    carbon_graph[idx2].add(idx1)
        
        # A terminal carbon (on the sugar backbone) will have only one connection.
        terminal_aldehyde_found = False
        # For each match of the aldehyde group, check if the carbon is terminal.
        for match in mol.GetSubstructMatches(aldehyde_pattern):
            # The match returns a tuple with the aldehyde carbon as the first atom.
            carbon_idx = match[0]
            if carbon_idx in carbon_graph and len(carbon_graph[carbon_idx]) == 1:
                terminal_aldehyde_found = True
                break
                
        if not terminal_aldehyde_found:
            return False, "Acyclic structure: aldehyde group is not located on a terminal carbon of the sugar backbone."
        
        return True, "Open-chain aldopentose: 5-carbon backbone with terminal aldehyde group confirmed."
    
    else:
        # For cyclic forms, even though the explicit aldehyde group might be absent,
        # the structure is considered an aldopentose if it is cyclic and does not contain lactone functionality.
        if mol.HasSubstructMatch(lactone_pattern):
            return False, "Cyclic structure contains lactone (ester) functionality; not classified as an aldopentose."
        return True, "Cyclized aldopentose: cyclic structure with correct 5-carbon backbone and no disqualifying lactone functionality."

# (Optional testing calls)
if __name__ == "__main__":
    # Example aldopentose SMILES (aldehydo-L-xylose)
    test_smiles = "[H]C(=O)[C@@H](O)[C@H](O)[C@@H](O)CO"
    result, reason = is_aldopentose(test_smiles)
    print(f"SMILES: {test_smiles}\nClassification: {result}\nReason: {reason}")