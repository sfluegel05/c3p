"""
Classifies: CHEBI:46722 carbonate ester
"""
#!/usr/bin/env python3
"""
Classifies: Carbonate ester

Defined as: Any derivative of carbonic acid (O=C(OH)2) in which one or both of the hydroxyl hydrogens 
have been replaced by organyl (carbon‐containing) groups. This includes both acyclic and cyclic carbonate 
esters (including monoesters where at least one –OH has been replaced), but not carbonic acid itself.

The strategy implemented here is:
  1. Loop over all carbon atoms that are only bonded to oxygen atoms and have exactly 3 such neighbors.
  2. Among these, require that exactly one bond is a double bond (the carbonyl) with an oxygen 
     that is “simple” (i.e. it has no extra heavy neighbours).
  3. For the one or two single‐bonded oxygens, require that at least one is “substituted” – that is, 
     it is connected to an additional carbon atom which is “organic” (we use a heuristic: the carbon is sp3 
     or has an explicit hydrogen).
     
If a candidate is found the molecule is classified as a carbonate ester (or monoester) moiety.
Note:
  This heuristic approach may still miss very “exotic” cases or may flag some borderline cases.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule contains a carbonate ester (or carbonate monoester) moiety based on its SMILES.
    The function examines each carbon atom and checks whether it is connected exclusively to oxygen atoms,
    has exactly three such neighbors (one double-bond and two single bonds) and that at least one of the 
    single-bonded oxygen atoms is substituted by an organyl group (i.e. attached to a carbon that is aliphatic).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a carbonate ester moiety is found, False otherwise.
        str: Explanation of the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Loop over candidate central carbons.
    for atom in mol.GetAtoms():
        # Candidate carbon must be atomic number 6.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Collect all neighbors that are oxygen.
        oxygen_neighbors = []
        oxygen_bonds = []
        valid_candidate = True
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                oxygen_neighbors.append(nbr)
                oxygen_bonds.append(bond)
            else:
                # If the candidate carbon is attached to a non-oxygen atom,
                # it is unlikely to be part of an isolated carbonate group.
                valid_candidate = False
                break
        if not valid_candidate:
            continue

        # For a carbonic acid derivative, expect exactly 3 oxygen neighbours.
        if len(oxygen_neighbors) != 3:
            continue
        
        # Identify the double-bonded oxygen (the carbonyl) and the single-bonded oxygens.
        dbl_bond_count = 0
        carbonyl_ok = False
        single_oxygens = []
        for bond in oxygen_bonds:
            o_atom = bond.GetOtherAtom(atom)
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                dbl_bond_count += 1
                # Check that the double-bonded oxygen is simple. In many cases we expect it not to be further bonded.
                # (Allow a degree of 1 or 2 if the extra connection is to a hydrogen.)
                if o_atom.GetDegree() > 2:
                    carbonyl_ok = False
                    break
                else:
                    carbonyl_ok = True
            elif bond.GetBondType() == rdchem.BondType.SINGLE:
                single_oxygens.append(o_atom)
            else:
                # If bond type is not single/double, skip candidate.
                carbonyl_ok = False
                break
        if dbl_bond_count != 1 or not carbonyl_ok:
            continue

        # For the carbonate moiety, at least one single-bonded oxygen must be substituted by an organyl group.
        # We define "organyl substitution" here as: the oxygen (other than the one bound to the candidate carbon) 
        # is attached to at least one carbon that appears to be aliphatic (sp3 or bearing explicit hydrogens).
        substitution_found = False
        for o_atom in single_oxygens:
            # There may be cases where the oxygen is part of a ring; we rely on its explicit bonds.
            for nbr in o_atom.GetNeighbors():
                # Exclude back-edge to candidate carbon.
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    # Check if the attached carbon is aliphatic: either SP3 or it carries explicit hydrogens.
                    if (nbr.GetHybridization() == rdchem.HybridizationType.SP3) or (nbr.GetTotalNumHs() > 0):
                        substitution_found = True
                        break
            if substitution_found:
                break
        
        if substitution_found:
            return True, f"Found carbonate ester moiety at carbon atom with index {atom.GetIdx()}."
    
    return False, "No carbonate ester moiety with necessary organyl substitution found."

# Example testing block (can be removed when using as a module)
if __name__ == "__main__":
    # A set of examples that were used in previous testing.
    examples = [
        ("4,5-Dichloro-1,3-dioxolan-2-one", "ClC1OC(OC1Cl)=O"),
        ("(S)-tolpyralate", "C=1C(=C(C(=C(C1)C(=O)C=2C=NN(C2O[C@@H](OC(OC)=O)C)CC)C)OCCOC)S(C)(=O)=O"),
        ("Ethylene carbonate", "O1CCOC1=O"),
        ("Ethyl 3-methylbut-3-enyl carbonate", "O(CCC(C)=C)C(OCC)=O"),
        ("Stagonosporyne D", "O=C1O[C@H]2[C@@](O)(C#CC(=C)C)[C@@H](O)CC[C@H]2O1"),
        ("spiropidion", "CCOC(=O)OC1=C(C(=O)N(C)C11CCN(CC1)OC)C1=C(C)C=C(Cl)C=C1C"),
        ("fluoren-9-ylmethyl hydrogen carbonate", "[H]C1(COC(O)=O)c2ccccc2-c2ccccc12"),
        ("Diethyl carbonate", "O(CC)C(OCC)=O"),
        ("(Z)-hex-3-en-1-yl methyl carbonate", "CC/C=C\\CCOC(OC)=O"),
        ("Stagonosporyne E", "O=C1O[C@@H]2C(=C(CO)[C@H]([C@@H]2O1)(O)C#CC(=C)C)O"),
        ("5-(undeca-3,5,7-trien-1-yl)-2-furyl hydrogen carbonate", "CCCC(=C(C(=C(C(=C(CCC1=CC=C(O1)OC(=O)O)[H])[H])[H])[H])[H])[H]"),
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