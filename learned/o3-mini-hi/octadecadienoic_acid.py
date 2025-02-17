"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid – a straight-chain C18 polyunsaturated fatty acid having exactly 2 carbon-carbon double bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is defined as any straight-chain C18 polyunsaturated fatty acid with exactly two 
    non-aromatic carbon-carbon double bonds, and must have a terminal carboxylic acid group.
    
    The algorithm works as follows:
    1. Parse the SMILES. If invalid, return failure.
    2. Ensure the molecule has a carboxylic acid functionality.
    3. Verify that the molecule has a total of 18 carbon atoms.
    4. Identify the carboxyl carbon from the acid group and then attempt to "walk" along a connected carbon chain.
       • From the acid carbon, follow the unique carbon neighbor (if available) to get the main chain.
       • At each step, ensure that there is no branching: each atom (except the chain endpoints) must have exactly 
         two carbon neighbors that are in the chain.
    5. Confirm that the total length of the chain is 18 atoms.
    6. Count the double bonds along the chain (the connectivity between successive chain atoms); 
       exactly 2 non-aromatic C=C bonds must be present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as octadecadienoic acid, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Look for a carboxylic acid substructure.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Missing carboxylic acid functionality"
    # Use the first acid match. In the SMARTS, the first atom is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Step 2: Count the total number of carbon atoms in the molecule.
    all_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(all_carbons) != 18:
        return False, f"Expected 18 carbon atoms but found {len(all_carbons)}"
    
    # Step 3: Identify the main chain starting from the acid carbon.
    # The acid carbon should have exactly one carbon neighbor that continues the fatty acid chain.
    chain = []
    chain_set = set()  # To keep track of carbon atoms in the discovered chain
    chain.append(acid_carbon)
    chain_set.add(acid_carbon.GetIdx())
    
    # Find the unique carbon neighbor of the acid carbon (excluding oxygens)
    acid_neighbors = [n for n in acid_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(acid_neighbors) != 1:
        return False, "Carboxyl carbon does not have exactly one carbon neighbor; not a straight-chain acid"
    current = acid_neighbors[0]
    chain.append(current)
    chain_set.add(current.GetIdx())
    prev = acid_carbon
    
    # Walk along the chain
    while True:
        # Get carbon neighbors of current, excluding the atom we came from.
        neighbors = [n for n in current.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != prev.GetIdx()]
        if len(neighbors) == 0:
            break  # reached the end of chain
        if len(neighbors) > 1:
            return False, f"Carbon atom with index {current.GetIdx()} shows branching in the carbon chain"
        # Move to the next neighbor
        next_atom = neighbors[0]
        # If we already visited this atom, then cycle (should not happen in a straight chain)
        if next_atom.GetIdx() in chain_set:
            break
        chain.append(next_atom)
        chain_set.add(next_atom.GetIdx())
        prev, current = current, next_atom

    if len(chain) != 18:
        return False, f"Main carbon chain length is {len(chain)} instead of 18; indicates branching or missing carbons"
    
    # Verify linear connectivity: endpoints should have one carbon neighbor within the chain; internal atoms exactly 2.
    for i, atom in enumerate(chain):
        # Count how many neighbors of this atom are also in the chain.
        cnt = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in chain_set:
                cnt += 1
        if i == 0 or i == len(chain)-1:
            if cnt != 1:
                return False, f"Chain endpoint atom index {atom.GetIdx()} has {cnt} neighbors within chain; should have 1"
        else:
            if cnt != 2:
                return False, f"Internal chain atom index {atom.GetIdx()} has {cnt} neighbors within chain; should have 2"
    
    # Step 4: Count the C=C bonds along the chain. They are the bonds connecting consecutive atoms in the chain.
    double_bonds = 0
    for i in range(len(chain)-1):
        bond = mol.GetBondBetweenAtoms(chain[i].GetIdx(), chain[i+1].GetIdx())
        if bond is None:
            return False, f"Missing bond between chain atoms at indices {chain[i].GetIdx()} and {chain[i+1].GetIdx()}"
        # Count bond as double if it is a double bond and non-aromatic.
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and (not bond.GetIsAromatic()):
            double_bonds += 1
    if double_bonds != 2:
        return False, f"Found {double_bonds} carbon-carbon double bonds along the main chain; exactly 2 required"
    
    return True, "Molecule is a straight-chain C18 fatty acid with exactly 2 C=C double bonds and a carboxylic acid group"

# Example testing (run only if this file is executed as a script):
if __name__ == "__main__":
    test_smiles = [
        "CCCCCC\\C=C/C=C\\CCCCCCCC(=O)O",  # 9Z,11Z-octadecadienoic acid (correct)
        "C(=C\\C/C=C\\CCCCCO)\\CCCCCCCC(=O)O",  # 18-hydroxylinoleic acid (should fail due to branching)
        "OC(=O)CCCCCCC\\C=C/C=C\\CCCCCC"  # 9Z,11Z-octadecadienoic acid alternative
    ]
    for sm in test_smiles:
        result, reason = is_octadecadienoic_acid(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")