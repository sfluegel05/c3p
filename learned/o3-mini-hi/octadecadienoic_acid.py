"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid – a straight-chain C18 polyunsaturated fatty acid 
having exactly 2 non-aromatic C=C double bonds and a terminal carboxylic acid group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) is an octadecadienoic acid.
    An octadecadienoic acid is defined as a straight-chain fatty acid that:
      • Contains exactly 18 carbon atoms,
      • Possesses a terminal carboxylic acid group ([CX3](=O)[OX2H1]),
      • Has a linear (unbranched) chain connected from the acid carbon (the carbonyl in COOH) 
        such that no chain atom bears an extra carbon substituent, and 
      • Contains exactly 2 non-aromatic carbon-carbon double bonds along the consecutive carbon chain.
    
    The algorithm is as follows:
    1. Parse the SMILES. If invalid, fail.
    2. Look for a carboxylic acid substructure using SMARTS "[CX3](=O)[OX2H1]". If not found, fail.
    3. Verify that the total number of carbons in the molecule is exactly 18.
    4. Starting with the carboxyl carbon (the first atom in the SMARTS match), determine its unique 
       carbon neighbor (it should have exactly one – otherwise the chain is not linear).
    5. “Walk” the chain by always taking the unique carbon neighbor (excluding the one we came from)
       until no further carbon neighbor is found. (This “path” is assumed to be the main chain.)
    6. Check that the length of the obtained chain is exactly 18. (This excludes any molecules with 
       extra carbons (branched or cyclic) on the main chain.)
    7. For each carbon atom in the chain, verify that *all* its carbon neighbors are in the chain.
       (This step detects branch points.)
    8. Examine the bonds along consecutive chain atoms and count the number of double bonds that are 
       non-aromatic; exactly 2 are required.
       
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple. The boolean is True if the molecule is an octadecadienoic acid 
                     and False otherwise, while str explains the reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check for carboxylic acid functionality.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Missing carboxylic acid functionality"
    
    # Use first match: by convention, the acid carbon (the carbonyl carbon) is the first atom.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Step 2: Verify that the total number of carbons is exactly 18.
    all_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(all_carbons) != 18:
        return False, f"Expected 18 carbon atoms but found {len(all_carbons)}"
    
    # Step 3: Build the main chain starting from the acid carbon.
    # The acid carbon should have exactly one carbon neighbor that continues the fatty acid chain.
    acid_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(acid_neighbors) != 1:
        return False, "Carboxyl carbon does not have exactly one carbon neighbor; not a straight-chain acid"
    
    chain = [acid_carbon]  # list of atoms in the main chain
    chain_indices = {acid_carbon.GetIdx()}
    
    # Add the unique neighbor of acid carbon.
    next_atom = acid_neighbors[0]
    chain.append(next_atom)
    chain_indices.add(next_atom.GetIdx())
    
    prev_atom = acid_carbon
    current_atom = next_atom

    # Walk along the chain:
    while True:
        # Get all carbon neighbors (atomic number == 6) not equal to prev_atom.
        nbr_carbons = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx()]
        if len(nbr_carbons) == 0:
            break  # reached the chain end
        if len(nbr_carbons) > 1:
            # This carbon has more than one carbon neighbor (not counting the one we came from)
            # which indicates a branching in the backbone.
            return False, f"Carbon atom with index {current_atom.GetIdx()} shows branching in the carbon chain"
        # There is exactly one possible next atom.
        next_atom = nbr_carbons[0]
        if next_atom.GetIdx() in chain_indices:
            # Already visited this carbon; would lead to a cycle (should not happen in a straight chain).
            break
        chain.append(next_atom)
        chain_indices.add(next_atom.GetIdx())
        prev_atom, current_atom = current_atom, next_atom

    # Step 4: Verify that the chain length is exactly 18.
    if len(chain) != 18:
        return False, f"Main carbon chain length is {len(chain)} instead of 18; indicates branching, cycle, or extra/missing carbons"

    # Step 5: Ensure that the chain is “pure” – no chain atom has any extra carbon neighbors that are not in the chain.
    for atom in chain:
        # List all carbon neighbors of the atom.
        carbon_nbrs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        for nbr in carbon_nbrs:
            if nbr.GetIdx() not in chain_indices:
                return False, f"Carbon atom with index {atom.GetIdx()} has an extra carbon substituent (branching detected)"

    # Step 6: Count double bonds (non-aromatic) between consecutive atoms in the chain.
    double_bond_count = 0
    for i in range(len(chain)-1):
        bond = mol.GetBondBetweenAtoms(chain[i].GetIdx(), chain[i+1].GetIdx())
        if bond is None:
            return False, f"Missing bond between chain atoms at indices {chain[i].GetIdx()} and {chain[i+1].GetIdx()}"
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and (not bond.GetIsAromatic()):
            double_bond_count += 1
    if double_bond_count != 2:
        return False, f"Found {double_bond_count} non-aromatic C=C bonds along the main chain; exactly 2 are required"

    return True, "Molecule is a straight-chain C18 fatty acid with exactly 2 C=C double bonds and a carboxylic acid group"

# For testing purposes (only executed if the module is run as a script)
if __name__ == "__main__":
    # Some test SMILES strings (these include examples from the successful and failing cases).
    test_smiles = [
        # True positives:
        "CCCCCC\\C=C/C=C/[C@H](O)CCCCCCCC(O)=O",  # 9(R)-HODE
        "OC(=O)CCCCCCC\\C=C/C=C\\CCCCCC",           # 9Z,11Z-octadecadienoic acid
        # False positives (should be rejected for branching):
        "C(=C\\C/C=C\\CCCCCO)CCCCCCCC(=O)O",         # 18-hydroxylinoleic acid
        "OC(=O)CCCCC/C=C/C=C\\CCCCCCCC",             # 7-trans,9-cis-octadecadienoic acid
        # False negative example (should reject because total carbons != 18):
        "O(C(CCCCCCC(O)=O)/C=C/C(=O)CCCCCCCCC(O)=O",  # (11E)-13-hydroxy-10-oxo-11-octadecenoic acid (has extra carbon)
    ]
    for sm in test_smiles:
        result, reason = is_octadecadienoic_acid(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")