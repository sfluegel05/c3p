"""
Classifies: CHEBI:46722 carbonate ester
"""
#!/usr/bin/env python3
"""
Classifies: Carbonate ester

Defined as: Any derivative of carbonic acid (O=C(OH)2) in which one or both of the hydroxyl 
hydrogens have been replaced by organyl (carbon‐containing) groups. 
That is, the functional group appears as –OC(=O)OR or –OC(=O)(OH) (when only one H is replaced)
or in cyclic variants.
 
The strategy is to look for a candidate carbon that is connected to exactly three oxygen atoms:
  • exactly one of these oxygens is connected via a double bond (the carbonyl oxygen), and 
  • the other two are connected via single bonds.
Then for at least one of the single‐bonded oxygens, we require that the oxygen is “substituted” by 
an organic group (i.e. its other neighbor is carbon rather than just hydrogen). This avoids classifying 
carbonic acid (where both –OH groups remain) as a carbonate ester.
"""

from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule contains a carbonate ester (or carbonate monoester) moiety based on its SMILES.
    The function scans through each carbon atom in the molecule and checks if it has exactly three attached
    oxygen atoms – one via a double bond and two via single bonds. Then it further checks that at least one 
    of the single-bonded oxygen atoms is substituted with an organyl group (i.e., attached to a carbon atom)
    rather than being only –OH (or similar).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a carbonate ester moiety is found, False otherwise.
        str: Explanation of the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Loop over all carbon atoms to identify a carbonate carbon candidate.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # Only carbon atoms can be the central carbon of a carbonate group.
        
        # Get all neighboring oxygens attached to this carbon.
        oxygen_neighbors = []
        oxygen_bonds = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                oxygen_neighbors.append(nbr)
                oxygen_bonds.append(bond)
        
        # In a carbonate moiety (CO3), the carbon must be attached to exactly three oxygen atoms.
        if len(oxygen_neighbors) != 3:
            continue
        
        # Count the number of double bonds to oxygen and collect the single bonds.
        dbl_count = 0
        single_oxygens = []
        for bond in oxygen_bonds:
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                dbl_count += 1
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                single_oxygens.append(bond.GetOtherAtom(atom))
        
        # The candidate must have exactly one double-bonded oxygen (the carbonyl) and two single-bonded oxygens.
        if dbl_count != 1 or len(single_oxygens) != 2:
            continue
        
        # For each single-bonded oxygen, check for substitution.
        # We require that at least one of these oxygens is attached (other than to the carbonate carbon)
        # to a carbon atom (organyl group). This ensures we are not simply looking at –OH.
        for o_atom in single_oxygens:
            # Get neighbors of oxygen excluding our candidate carbon.
            other_neighbors = [n for n in o_atom.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
            # Check if any of the other neighbors is carbon (atomic number 6)
            if any(n.GetAtomicNum() == 6 for n in other_neighbors):
                return True, f"Found carbonate ester moiety at carbon atom with index {atom.GetIdx()}."
        
        # If neither oxygen is substituted, then the moiety is likely just carbonic acid.
    
    return False, "No carbonate ester moiety with necessary organyl substitution found."


# Example testing of the function with a few representative SMILES.
if __name__ == "__main__":
    examples = [
        ("4,5-Dichloro-1,3-dioxolan-2-one", "ClC1OC(OC1Cl)=O"),
        ("(S)-tolpyralate", "C=1C(=C(C(=C(C1)C(=O)C=2C=NN(C2O[C@@H](OC(OC)=O)C)CC)C)OCCOC)S(C)(=O)=O"),
        ("Ethylene carbonate", "O1CCOC1=O"),
        ("Ethyl 3-methylbut-3-enyl carbonate", "O(CCC(C)=C)C(OCC)=O"),
        ("Stagonosporyne D", "O=C1O[C@H]2[C@@](O)(C#CC(=C)C)[C@@H](O)CC[C@H]2O1"),
        ("spiropidion", "CCOC(=O)OC1=C(C(=O)N(C)C11CCN(CC1)OC)C1=C(C)C=C(Cl)C=C1C"),
        ("Diethyl carbonate", "O(CC)C(OCC)=O"),
        ("(Z)-hex-3-en-1-yl methyl carbonate", "CC/C=C\\CCOC(OC)=O"),
        ("Stagonosporyne E", "O=C1O[C@@H]2C(=C(CO)[C@H]([C@@H]2O1)(O)C#CC(=C)C)O"),
        ("diphenyl carbonate", "O=C(Oc1ccccc1)Oc1ccccc1"),
        ("(3aR,4S,5R,7aS)-4,5-dihydroxy-6-methyl-3a,4,5,7a-tetrahydro-1,3-benzodioxol-2-one", "O=C1O[C@H]2C=C(C)[C@H]([C@@H]([C@H]2O1)O)O"),
        ("methyl 4-nitrobenzyl carbonate", "C=1(C=CC(=CC1)[N+]([O-])=O)COC(OC)=O"),
        ("Phomoxin C", "O=C1O[C@H]2C(=C(CO)[C@H]([C@@H]([C@H]2O1)O)O)/C=C/CCCCC"),
        ("Phomoxin", "O=C1O[C@H]2C(=C(CO)[C@@H]([C@@H]([C@H]2O1)O)O)/C=C/CCCCC"),
        ("monoethyl carbonate", "CCOC(O)=O"),
        ("dimethyl carbonate", "COC(=O)OC"),
        ("(3aS,4R,5S,7aR)-4,5-dihydroxy-7-methyl-3a,4,5,7a-tetrahydrobenzo[1,3]dioxol-2-one", "O=C1O[C@@H]2C(=C[C@@H]([C@H]([C@@H]2O1)O)O)C")
    ]
    
    for name, smi in examples:
        result, reason = is_carbonate_ester(smi)
        print(f"{name}: {result} ({reason})")