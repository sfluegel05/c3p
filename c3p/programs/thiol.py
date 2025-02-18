"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol – An organosulfur compound in which a thiol group (-SH or its equivalent, e.g. a thioether substitution)
is attached to a carbon atom of any aliphatic or aromatic moiety.
However, when the S-containing group is part of an amino acid or peptide backbone (i.e. the carbon attached to S resembles an 
α–carbon with both an amine and a carbonyl-bearing neighbor), the molecule should not be classified as a simple thiol.
"""

from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a (nonpeptidic) thiol based on its SMILES string.
    
    The strategy is:
      1. Parse the molecule and add explicit hydrogens.
      2. Loop over all sulfur atoms that are not "oxidized" (i.e. not double bonded to oxygen).
      3. For each S atom attached to at least one carbon, check each S–C bond:
           - We use a heuristic to detect an amino acid α–carbon environment.
             Specifically, if the carbon attached to S (call it C_thiol) has a neighboring nitrogen (N) and also a neighboring
             carbon which is double-bonded to an oxygen (i.e. a carbonyl group), then we consider this bond to be peptide‐like.
      4. If at least one S–C bond is found that is NOT in a peptide-like environment, we return True.
      5. Otherwise, (either no S–C bonds found or all are peptide-like) we return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if at least one S atom is attached to a carbon outside an amino acid/peptide-like environment.
        str: Explanation / Reason for the classification.
    """
    # Step 1: Parse molecule and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Flag to see if at least one candidate S–C bond is found.
    candidate_found = False
    
    # Loop over all atoms: we are looking for sulfur atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # skip non-sulfur atoms
        # Step 2: Skip sulfur if it is "oxidized" (i.e. has a double bond to an oxygen)
        oxidized = any(bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(atom).GetAtomicNum() == 8 
                       for bond in atom.GetBonds())
        if oxidized:
            continue  # skip oxidized sulfur groups (sulfoxides, sulfones, etc.)
        
        # Now examine all neighbors of this sulfur.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue  # we only care about S attached to carbon
            # At this stage we have a candidate S-C bond.
            # Now check if the carbon (which we call C_thiol) appears to be part of a peptide-like environment.
            peptide_like = False
            
            # Heuristic: if the C_thiol has at least one neighbor that is nitrogen (N)
            # and at least one other neighbor that is a carbon (C) which is double-bonded to an oxygen.
            has_N = any(n.GetAtomicNum() == 7 for n in nbr.GetNeighbors() if n.GetIdx() != atom.GetIdx())
            has_carbonyl = False
            for second in nbr.GetNeighbors():
                if second.GetIdx() == atom.GetIdx():
                    continue
                if second.GetAtomicNum() == 6:
                    # Examine bonds of this neighboring carbon to see if it has a double bond to oxygen.
                    for bond in second.GetBonds():
                        if bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(second).GetAtomicNum() == 8:
                            has_carbonyl = True
                            break
                if has_carbonyl:
                    break
            
            if has_N and has_carbonyl:
                peptide_like = True
            
            # If this S–C bond is not in a peptide-like environment, we consider the molecule to contain a simple thiol.
            if not peptide_like:
                return True, "Molecule contains a sulfur group (e.g. -SH or thioether substitution) attached to a carbon atom outside a peptide backbone"
    # If we finish the loop without returning, then either no S–C bond was found or every candidate was peptide-like.
    return False, "All detected S groups are associated with peptide backbones or no valid thiol-like group (S attached to C) found"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # A few examples to illustrate behavior:
    test_molecules = {
        "2-methoxy-4-(methylthio)-N-pyridin-4-ylbenzamide": "COC1=C(C=CC(=C1)SC)C(=O)NC2=CC=NC=C2",
        "(S)-fluoxytioconazole (thiol)": "O[C@](CN1N=CN=C1S)(C1=CC=C(F)C=C1F)C(F)(F)C1=CC=C(OC2=CC=C(C=C2)C#N)C=N1",
        "coenzyme A": "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        "homocysteine": "NC(CCS)C(O)=O",  # should be rejected as peptide-like
        "mercaptoethanol": "OCCS",  # simple thiol
    }
    
    for name, smi in test_molecules.items():
        result, reason = is_thiol(smi)
        print(f"{name}:\n  SMILES: {smi}\n  Classified as thiol? {result}\n  Reason: {reason}\n")