"""
Classifies: CHEBI:46722 carbonate ester
"""
#!/usr/bin/env python3
"""
Classifies: Carbonate ester
Defined as: Any carbonate (derived from carbonic acid, O=C(OH)2)
in which the hydroxyl hydrogens have been replaced by organyl groups.
That is, the functional group should appear as –OC(=O)OR (or cyclic variants)
with both -OH replaced.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES.
    The function identifies a carbon atom that is bonded to a double bonded oxygen
    (C=O) and two single-bonded oxygen atoms. Each of those oxygen atoms (the ester oxygens)
    must be substituted with an organic (non-hydrogen) group rather than a proton.
    
    Note: Some provided examples (like monoethyl carbonate) are only partially esterified.
    Here we require that both –OH groups are replaced as per the definition.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a carbonate ester moiety is found, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Iterate over all atoms looking for a potential carbonate carbon.
    # We expect a carbonate carbon to be sp2 hybridized with three oxygen neighbors:
    # one via a double bond (C=O) and two via single bonds (C–O).
    for atom in mol.GetAtoms():
        # Focus on carbon atoms only.
        if atom.GetAtomicNum() != 6:
            continue
        
        doubleO = None
        singleOs = []
        
        # Examine bonds around the carbon atom.
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            # We are interested in bonds to oxygen only.
            if nbr.GetAtomicNum() != 8:
                continue
            
            # Check bond type.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Expect only one C=O.
                if doubleO is None:
                    doubleO = nbr
                else:
                    # More than one double bonded oxygen? Not typical for a carbonate.
                    doubleO = None
                    break
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                singleOs.append(nbr)
        
        # We need exactly one double-bonded oxygen and exactly two single-bonded oxygen atoms.
        if doubleO is None or len(singleOs) != 2:
            continue  # not a candidate carbonate group
        
        # Check each of the single-bonded oxygens to ensure they are substituted.
        # They should be not just –OH; i.e., besides the carbonyl carbon, they must be attached to
        # at least one non-hydrogen atom.
        all_substituted = True
        for o_atom in singleOs:
            # Get neighbors of this oxygen, excluding the carbonyl carbon.
            neighbors = [n for n in o_atom.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
            if not neighbors:
                all_substituted = False
                break
            # For a proper ester oxygen in a diester there should be exactly one other neighbor.
            # Check that this neighbor is not a hydrogen.
            # (If the neighbor is hydrogen, then it is a hydroxyl group.)
            found_organyl = any(n.GetAtomicNum() != 1 for n in neighbors)
            if not found_organyl:
                all_substituted = False
                break
        
        if all_substituted:
            return True, f"Found carbonate ester moiety at carbon atom with index {atom.GetIdx()}."
    
    return False, "No properly substituted carbonate (–OC(=O)OR) moiety found that matches the definition."

# You can test the function with some SMILES examples:
if __name__ == "__main__":
    examples = [
        ("4,5-Dichloro-1,3-dioxolan-2-one", "ClC1OC(OC1Cl)=O"),
        ("Ethylene carbonate", "O1CCOC1=O"),
        ("Diethyl carbonate", "COC(=O)OC"),
        ("Monoethyl carbonate", "CCOC(O)=O"),
        ("diphenyl carbonate", "O=C(Oc1ccccc1)Oc1ccccc1")
    ]
    
    for name, smi in examples:
        result, reason = is_carbonate_ester(smi)
        print(f"{name}: {result} ({reason})")