"""
Classifies: CHEBI:46722 carbonate ester
"""
#!/usr/bin/env python3
"""
Classifies: Carbonate ester

Defined as: Any derivative of carbonic acid (O=C(OH)2) in which at least one of the hydroxyl hydrogens 
has been replaced by an organyl group. That is, the functional group appears as –OC(=O)OR 
or cyclic variants. (Note: Fully esterified carbonates have both –OH replaced. Here we allow partial substitution 
to catch e.g. monoesters as well.)

The strategy is to look for a candidate carbon that is connected to exactly three oxygen atoms:
  • one via a double bond (C=O) and 
  • two via single bonds (C–O).
Then we check that these two “ester” oxygens are not both only bound to a hydrogen (in which case the moiety would be carbonic acid).
This approach distinguishes many carbonate esters from other moieties.
"""

from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule contains a carbonate ester (or carbonate monoester) moiety based on its SMILES.
    The function scans through all carbon atoms and checks whether one has exactly one double-bonded oxygen 
    and two single-bonded oxygen atoms. If both single-bonded oxygens are unsubstituted (i.e. only have a hydrogen attached 
    besides the carbonyl carbon), then the moiety is just carbonic acid and is not classified as an ester.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a carbonate ester moiety is found, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Iterate over all atoms to find a candidate carbonate-carbon.
    # Ideally, the candidate carbon has three neighbors (degree 3): one double bonded oxygen and two single bonded oxygens.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # To be a carbonate carbon the connectivity should be exactly three.
        if atom.GetDegree() != 3:
            continue

        double_bonded_oxygens = []
        single_bonded_oxygens = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() != 8:
                continue  # only care about bonds to oxygen
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                double_bonded_oxygens.append(nbr)
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                single_bonded_oxygens.append(nbr)
        
        # Check required connectivity: one double-bond and exactly two single bonds to oxygen.
        if len(double_bonded_oxygens) != 1 or len(single_bonded_oxygens) != 2:
            continue
        
        # For each single-bonded oxygen, check its other connections (besides the candidate carbon).
        # If both oxygens are only bound to hydrogen (i.e. –OH) then this is just carbonic acid.
        substitution_flags = []
        for o_atom in single_bonded_oxygens:
            other_neighbors = [n for n in o_atom.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
            # Mark the oxygen as substituted if at least one neighbor is not hydrogen.
            is_substituted = any(n.GetAtomicNum() != 1 for n in other_neighbors)
            substitution_flags.append(is_substituted)
        # If at least one oxygen is substituted by an organyl group, we classify it as a carbonate ester.
        if substitution_flags[0] or substitution_flags[1]:
            return True, f"Found carbonate ester moiety at carbon atom with index {atom.GetIdx()}."
    
    return False, "No properly substituted carbonate (–OC(=O)OR or –OC(=O)(OH)) moiety found that matches the definition."

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
        ("(3aS,4R,5S,7aR)-4,5-dihydroxy-7-methyl-3a,4,5,7a-tetrahydrobenzo[1,3]dioxol-2-one", "O=C1O[C@@H]2C(=C[C@@H]([C@H]([C@@H]2O1)O)O)C"),
        # Additional examples can be added here.
    ]
    
    for name, smi in examples:
        result, reason = is_carbonate_ester(smi)
        print(f"{name}: {result} ({reason})")