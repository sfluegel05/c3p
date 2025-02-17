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
  1. Loop over all carbon atoms and select those that are sp2 and have exactly three neighbours 
     (as required by the trigonal planar structure of a carbonate group).
  2. Among the attached oxygens, require that exactly one is double‐bonded and that oxygen is “simple”
     (has degree 1) – identifying the carbonyl oxygen.
  3. For the one or two single‐bonded oxygens, require that at least one has a neighbour (other than the candidate
     carbon) that is a carbon (atomic number 6). This indicates that an organyl substitution has replaced an –OH.
     
If such a candidate is found the molecule is classified as containing a carbonate ester (or carbonate monoester)
moiety.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule contains a carbonate ester (or carbonate monoester) moiety based on its SMILES.
    The function examines each carbon atom and checks whether it has three oxygen neighbors with
    the expected bonding pattern (one oxygen double‐bonded having degree 1, and at least one single‐bonded oxygen
    that is substituted by a carbon rather than being simply –OH).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a carbonate ester moiety is found, False otherwise.
        str: Explanation of the outcome.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Loop over each atom and look for a candidate carbonate carbon.
    for atom in mol.GetAtoms():
        # Candidate central atom must be carbon.
        if atom.GetAtomicNum() != 6:
            continue
        
        # For a carbonate group (derived from carbonic acid) the carbon is sp2 and bonded to exactly three atoms.
        if atom.GetHybridization() != rdchem.HybridizationType.SP2:
            continue
        if atom.GetDegree() != 3:
            continue
        
        # Collect oxygen neighbors and the bonds connecting them.
        oxygen_neighbors = []
        oxygen_bonds = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                oxygen_neighbors.append(nbr)
                oxygen_bonds.append(bond)
        
        # To be a carbonic acid derivative, the candidate carbon should have exactly 3 oxygen neighbours.
        if len(oxygen_neighbors) != 3:
            continue
        
        # Identify the oxygen that is double-bonded (the carbonyl oxygen) and those with single bonds.
        dbl_bond_count = 0
        carbonyl_found = False
        single_oxygens = []
        for bond in oxygen_bonds:
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                dbl_bond_count += 1
                # Check that the double-bonded oxygen is “simple” (typically degree 1)
                o_atom = bond.GetOtherAtom(atom)
                # In resonance the degree might be slightly higher, but usually degree 1 is expected.
                if o_atom.GetDegree() != 1:
                    carbonyl_found = False
                    break
                else:
                    carbonyl_found = True
            elif bond.GetBondType() == rdchem.BondType.SINGLE:
                single_oxygens.append(bond.GetOtherAtom(atom))
            else:
                # Skip if bond type is not single or double.
                continue
        # We expect exactly one double bond to oxygen.
        if dbl_bond_count != 1 or not carbonyl_found:
            continue
        
        # For the carbonate moiety the remaining single-bonded oxygens are either –OH or –OR.
        # At least one must be substituted (i.e. attached to a carbon) so that we are not simply looking at carbonic acid.
        substitution_found = False
        for o_atom in single_oxygens:
            # Get neighbors of oxygen excluding the candidate carbon.
            for nbr in o_atom.GetNeighbors():
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

# Example testing block (can be removed when using as a module)
if __name__ == "__main__":
    # A list of (name, SMILES) pairs to illustrate the behavior.
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