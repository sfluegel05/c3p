"""
Classifies: CHEBI:33566 catechols
"""
#!/usr/bin/env python3
"""
Classifies: Catechols (Any compound containing an o-diphenol component)
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol (contains an o-diphenol component)
    based on its SMILES string. It finds an aromatic ring that has at least one pair
    of adjacent carbon atoms both bearing a hydroxyl (-OH) substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains an o-diphenol (catechol) moiety, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups are explicitly defined
    mol = Chem.AddHs(mol)
    
    # Create a dictionary to mark aromatic carbons that have an -OH substituent.
    # We iterate over all atoms, and for aromatic carbons, we check neighbors.
    hydroxylated_aromatic = set()
    for atom in mol.GetAtoms():
        # Check if atom is carbon and aromatic
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            # Check neighbors for an -OH group (oxygen that is bonded to at least one hydrogen)
            for neighbor in atom.GetNeighbors():
                # Check if neighbor is oxygen
                if neighbor.GetAtomicNum() == 8:
                    # Look at the neighbors of this oxygen to see if it is an -OH (one hydrogen attached)
                    has_hydrogen = any(nbr.GetAtomicNum() == 1 for nbr in neighbor.GetNeighbors())
                    if has_hydrogen:
                        hydroxylated_aromatic.add(atom.GetIdx())
                        # We only need one -OH per aromatic carbon for our purpose
                        break
                        
    # Now, check for any bond connecting two aromatic carbons such that both are hydroxylated.
    # This qualifies as an o-diphenol (catechol) substructure.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and 
            a1.GetIsAromatic() and a2.GetIsAromatic()):
            if a1.GetIdx() in hydroxylated_aromatic and a2.GetIdx() in hydroxylated_aromatic:
                return True, "Contains o-diphenol (catechol) moiety on an aromatic ring"
                
    return False, "No adjacent hydroxyl groups on an aromatic ring found"

# You can uncomment the following lines to test the function with some examples:
# examples = [
#     "CC(C)c1cccc(O)c1O",    # 3-isopropylcatechol, expected True
#     "CCCC",                 # Simple alkane, expected False
# ]
# for smi in examples:
#     result, reason = is_catechols(smi)
#     print(f"SMILES: {smi} -> {result}: {reason}")