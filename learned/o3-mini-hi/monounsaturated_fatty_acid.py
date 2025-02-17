"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated Fatty Acid

Definition: Any fatty acid with one double (but not triple) bond in the long, acyclic fatty acid chain,
with all remaining carbon-carbon bonds being single. The molecule must contain a terminal carboxylic acid
group and have a sufficient number of carbons (here, at least 5) to qualify as a fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    
    Steps:
      1. Parse the SMILES string.
      2. Reject cyclic molecules.
      3. Identify a carboxylic acid group (protonated or deprotonated).
      4. Verify that the acid carbon (the carbonyl carbon) is terminal (i.e. it has exactly one carbon neighbor).
      5. Require the molecule to have at least 5 carbon atoms.
      6. Examine every bond:
            - Immediately reject if a triple bond between carbons is found.
            - For each C–C double bond, ensure that it is internal by checking that neither atom is too close 
              (within one bond) to the acid carbon (using the distance matrix) and that both atoms have at least 
              2 neighboring carbons (ensuring that the unsaturation is not terminal).
      7. Only if exactly one such double bond exists is the molecule classified as a monounsaturated fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): Tuple with True and a reason string if classification is positive;
                     otherwise, False and a reason message.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject cyclic molecules
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule is cyclic, which is not typical for a fatty acid chain"
    
    # Identify a carboxylic acid group using a SMARTS pattern that covers protonated and deprotonated forms.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Molecule does not contain a carboxylic acid group, thus not a fatty acid"

    # For simplicity, take the first match and identify the acid carbon (the carbon in C(=O)[O;H,-]).
    acid_match = acid_matches[0]
    acid_idx = acid_match[0]
    acid_atom = mol.GetAtomWithIdx(acid_idx)
    
    # The carboxylic acid's carbon should be terminal:
    carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Acid group is not terminal; the acid carbon has an unexpected bonding pattern"
    
    # Require that the molecule has a sufficient number of carbons (at least 5)
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 5:
        return False, "Molecule does not have enough carbon atoms to be considered a fatty acid"
    
    # Prepare a distance matrix to assess how far atoms are from the acid carbon
    dist_matrix = rdmolops.GetDistanceMatrix(mol)

    # Now count the carbon–carbon unsaturations.
    # We will count only double bonds (ignoring triple bonds) that are internal.
    unsat_double_count = 0
    for bond in mol.GetBonds():
        # Immediately reject if a triple bond between carbons is present
        if bond.GetBondType() == Chem.BondType.TRIPLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                return False, "Molecule has a triple bond between carbons, disqualifying it as a monounsaturated fatty acid"
        
        # Process double bonds between carbons
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6:
                continue  # skip if either atom is not carbon
            
            # Check that neither atom is too close to the acid carbon (i.e. must be at least 2 bonds away)
            if dist_matrix[acid_idx, a1.GetIdx()] <= 1 or dist_matrix[acid_idx, a2.GetIdx()] <= 1:
                continue  # likely part of the acid group region
            
            # Check that the unsaturation is internal.
            # For this, require that each atom in the double bond has at least two neighboring carbons.
            a1_c_nb = [nbr for nbr in a1.GetNeighbors() if nbr.GetAtomicNum() == 6]
            a2_c_nb = [nbr for nbr in a2.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(a1_c_nb) < 2 or len(a2_c_nb) < 2:
                continue  # terminal unsaturation; skip
            
            unsat_double_count += 1

    if unsat_double_count != 1:
        return False, f"Fatty acid chain unsaturation count is {unsat_double_count} rather than the required 1"

    return True, ("Molecule contains a terminal carboxylic acid group, is acyclic, has a sufficient carbon chain "
                  "length, and exactly one internal C–C double bond (with no triple bonds), classifying it as a "
                  "monounsaturated fatty acid")

# Example usage (uncomment to test):
# if __name__ == "__main__":
#     test_smiles = "CCCCCCCC\\C=C/CCCCCCC(O)=O"  # cis-8-heptadecenoic acid
#     result, reason = is_monounsaturated_fatty_acid(test_smiles)
#     print(result, reason)