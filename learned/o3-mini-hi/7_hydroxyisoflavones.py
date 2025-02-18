"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python
"""
Classifies: 7-hydroxyisoflavones 
Definition: A hydroxyisoflavone compound having a hydroxy group at the 7-position.
The function is_7_hydroxyisoflavones takes a SMILES string, verifies that the 
molecule contains an isoflavone-like bicyclic core and that one of the positions 
on the A-ring (adjacent to the heterocyclic oxygen) carries a free hydroxyl (-OH) group. 
Examples in the literature include structures such as luteone, wighteone, 7-hydroxyisoflavone, etc.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    It does so in two steps:
      1. Checks for an isoflavone-like core. Here we use a simplified SMARTS that matches 
         a bicyclic system with a pyran oxygen and a carbonyl group.
      2. Searches the core for a candidate 7-OH group. We assume that in a 7-hydroxyisoflavone 
         the bridging oxygen (the heterocyclic oxygen in the pyran ring) is bonded to an aromatic
         carbon (from the A-ring) that carries an –OH substituent.
         
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as a 7-hydroxyisoflavone, False otherwise.
        str: A reason message for the classification.
    """
    # Parse SMILES to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a simplified SMARTS pattern for an isoflavone core.
    # This pattern looks for a fused bicyclic system containing:
    #   - a pyran ring with a bridging oxygen (the 'o' atom) 
    #   - a ketone (carbonyl) group in the second ring.
    # The pattern “c1coc2cc(=O)cc(c2c1)” is a heuristic query.
    isoflavone_core = Chem.MolFromSmarts("c1coc2cc(=O)cc(c2c1)")
    if isoflavone_core is None:
        return False, "Error in SMARTS definition"
    
    # Check whether the molecule contains the isoflavone core
    core_matches = mol.GetSubstructMatches(isoflavone_core)
    if not core_matches:
        return False, "Isoflavone core not found"
    
    # For each matching core, try to find a hydroxyl group attached to the aromatic 
    # carbon of the A-ring that is adjacent to the heterocyclic oxygen.
    found7OH = False
    for match in core_matches:
        # The SMARTS pattern "c1coc2cc(=O)cc(c2c1)" is defined so that:
        #   match[1] corresponds to the heterocyclic oxygen in the pyran ring.
        try:
            bridging_oxygen = mol.GetAtomWithIdx(match[1])
        except Exception:
            continue
        
        # Look at neighbors of this bridging oxygen.
        # One neighbor should be the aromatic carbon belonging to the A-ring.
        for nbr in bridging_oxygen.GetNeighbors():
            # We choose only aromatic carbon neighbors:
            if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
                # Now, check if this aromatic carbon carries a hydroxyl substituent.
                # (We check the neighbors of the carbon for an oxygen that is not a carbonyl)
                for sub in nbr.GetNeighbors():
                    # Skip if the neighbor is the bridging oxygen itself:
                    if sub.GetIdx() == bridging_oxygen.GetIdx():
                        continue
                    if sub.GetAtomicNum() == 8:
                        # Exclude carbonyl oxygens (which are double-bonded to C)
                        # and require that the oxygen has at least one attached hydrogen.
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), sub.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                            continue
                        if sub.GetTotalNumHs() >= 1:
                            found7OH = True
                            break
                if found7OH:
                    break
        if found7OH:
            break
        
    if found7OH:
        return True, "Isoflavone core detected with hydroxy group at the presumed 7-position"
    else:
        return False, "Isoflavone core detected, but no hydroxy group found at the 7-position"

# Example usage (you can remove or comment these out in production):
if __name__ == "__main__":
    # Test with the 7-hydroxyisoflavone example
    smiles_examples = [
        "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O",  # 7-hydroxyisoflavone
        "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O"   # calycosin (has a 7-OH on the A-ring)
    ]
    for s in smiles_examples:
        result, reason = is_7_hydroxyisoflavones(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")