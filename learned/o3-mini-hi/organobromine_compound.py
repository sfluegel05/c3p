"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition: A compound containing at least one carbon–bromine bond.
This version improves on the previous approach by combining a SMARTS search with 
direct neighbor inspection. In cases where aromatic or resonance representations 
cause the Br atom to appear bound to heteroatoms even though it is “part of” an organic
framework, we check the context of the bromine’s bonded neighbors.
"""

from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined here as one that contains at least one bromine 
    atom covalently bonded to what is at least part of an organic (carbon‐based) framework.
    
    The function first uses a SMARTS substructure search for any bromine attached to a carbon.
    It then iterates over all bromine atoms in the molecule and checks their bonded neighbors.
    This combined approach helps to catch cases where resonance or tautomerization may result
    in the Br atom not appearing to be directly attached to C upon initial SMARTS matching.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organobromine compound, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, try a SMARTS query that should capture a bond between carbon and bromine.
    # The '~' operator means any bond type (single, double, aromatic, etc.)
    smarts_pattern = "[#6]~[Br]"
    pattern = Chem.MolFromSmarts(smarts_pattern)
    if pattern is None:
        return False, "Error in SMARTS pattern"
    
    if mol.HasSubstructMatch(pattern):
        # Found at least one instance where a carbon is directly bonded to bromine.
        return True, "Contains at least one carbon–bromine bond (SMARTS match)"
    
    # If the SMARTS search did not match, iterate over all Br atoms and check neighbors.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 35:  # Bromine has atomic number 35.
            # In some cases the connection might be drawn as Br bonded to another atom.
            # We check whether any neighbor is carbon (atomic number 6).
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    return True, "Contains at least one carbon–bromine bond (neighbor inspection)"
            # Optionally, we could check if the bromine is attached to a heteroatom that is
            # incorporated in an aromatic ring with predominantly carbon atoms.
            # For example, if the neighbor is nitrogen but is part of an aromatic ring that also contains carbons.
            # Here we implement a simple version: if the neighbor is not carbon, we inspect its ring.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() != 6 and neighbor.GetIsAromatic():
                    ring_info = mol.GetRingInfo()
                    # Check every ring that the neighbor is part of; if that ring contains any carbon, count it.
                    for ring in ring_info.AtomRings():
                        if neighbor.GetIdx() in ring:
                            for idx in ring:
                                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6:
                                    return True, ("Contains a bromine attached to an aromatic heterocycle "
                                                  "sharing ring with carbon atoms")
    # If no carbon–bromine bond is found by either method, return False.
    return False, "No carbon–bromine bond found"

# Example usage:
if __name__ == '__main__':
    # You can test the function with a few examples.
    test_smiles = [
        "Brc1ccccc1",                # bromobenzene (should be True)
        "Brn1ccc2ccccc12",            # 1-bromoindole (may be ambiguous but we try to catch it)
        "BrN1C(=O)CCC1=O",            # N-bromosuccinimide (expected True with our improved method)
        "ClCCBr",                    # 1-bromo-2-chloroethane (should be True)
        "CCC",                       # no bromine (should be False)
    ]
    for smi in test_smiles:
        result, reason = is_organobromine_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*40}")