"""
Classifies: CHEBI:46722 carbonate ester
"""
#!/usr/bin/env python3
"""
Classifies: Carbonate ester

Defined as: Any derivative of carbonic acid (O=C(OH)2) in which one or both of the hydroxyl hydrogens 
have been replaced by organyl (carbon‐containing) groups. This covers both acyclic and cyclic carbonate ester
(or monoester) moieties, but not carbonic acid itself.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_carbonate_ester(smiles: str):
    """
    Determines if the molecule contains a carbonate ester (or carbonate monoester) moiety.
    The algorithm looks for a central carbon that:
      - Is sp2 hybridized (as expected in a carbonyl carbon)
      - Has exactly three bonds, all to oxygen atoms.
      - Has exactly one double-bond to an oxygen (the carbonyl oxygen) that is not heavily substituted.
      - Has two single-bonded oxygens, at least one of which is substituted with any carbon-containing (organyl) group.
      
    Args:
        smiles (str): A SMILES string for the molecule.
        
    Returns:
        bool: True if at least one carbonate ester moiety is found, False otherwise.
        str: Explanation of the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    for atom in mol.GetAtoms():
        # Look only at carbon atoms
        if atom.GetAtomicNum() != 6:
            continue
        
        # Enforce that the candidate carbon should be sp2 (as in a carbonyl carbon)
        if atom.GetHybridization() != rdchem.HybridizationType.SP2:
            continue
        
        # Collect oxygen neighbors and check that there are exactly 3 bonds
        oxygen_neighbors = []
        oxygen_bonds = []
        # Also ensure degree (number of bonds) is exactly 3.
        if atom.GetDegree() != 3:
            continue
        
        valid_candidate = True
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:  # oxygen
                oxygen_neighbors.append(nbr)
                oxygen_bonds.append(bond)
            else:
                valid_candidate = False
                break
        if not valid_candidate or len(oxygen_neighbors) != 3:
            continue
        
        # Among these bonds, identify the one double bond (C=O) and the two single bonds.
        dbl_bond_count = 0
        carbonyl_ok = False
        single_oxygens = []
        for bond in oxygen_bonds:
            o_atom = bond.GetOtherAtom(atom)
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                dbl_bond_count += 1
                # For a carbonyl oxygen we allow a little extra connectivity if it is only to hydrogen or so;
                # here we check that it is not connected to an abundance of heavy atoms.
                if o_atom.GetDegree() > 2:
                    carbonyl_ok = False
                    break
                else:
                    carbonyl_ok = True
            elif bond.GetBondType() == rdchem.BondType.SINGLE:
                single_oxygens.append(o_atom)
            else:
                # If encountering unusual bond types, skip candidate.
                carbonyl_ok = False
                break
        
        if dbl_bond_count != 1 or not carbonyl_ok:
            continue
        
        # For the carbonate ester moiety, at least one single-bonded oxygen must be substituted
        # by an organyl group (any carbon; we do not restrict by hybridization so that aromatic carbons count)
        substitution_found = False
        for o_atom in single_oxygens:
            for nbr in o_atom.GetNeighbors():
                # Do not consider back edge to the candidate carbon.
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    substitution_found = True
                    break
            if substitution_found:
                break
        
        if substitution_found:
            return True, f"Found carbonate ester moiety at carbon atom with index {atom.GetIdx()}."
    
    return False, "No carbonate ester moiety with necessary organyl substitution found."

# Example testing block (can be removed when using this as a module)
if __name__ == "__main__":
    test_examples = [
        ("4,5-Dichloro-1,3-dioxolan-2-one", "ClC1OC(OC1Cl)=O"),
        ("(S)-tolpyralate", "C=1C(=C(C(=C(C1)C(=O)C=2C=NN(C2O[C@@H](OC(OC)=O)C)CC)C)OCCOC)S(C)(=O)=O"),
        ("Ethylene carbonate", "O1CCOC1=O"),
        ("Ethyl 3-methylbut-3-enyl carbonate", "O(CCC(C)=C)C(OCC)=O"),
        ("Stagonosporyne D", "O=C1O[C@H]2[C@@](O)(C#CC(=C)C)[C@@H](O)CC[C@H]2O1"),
        ("diphenyl carbonate", "O=C(Oc1ccccc1)Oc1ccccc1"),  # should now be accepted
        ("monoethyl carbonate", "CCOC(O)=O"),
        ("dimethyl carbonate", "COC(=O)OC")
    ]
    
    for name, smi in test_examples:
        res, reason = is_carbonate_ester(smi)
        print(f"{name}: {res} ({reason})")