"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen Compound
A compound containing at least one carbon–halogen bond (where X is a halogen atom).
Here we improve the previous approach by:
  1. First checking for a direct C–X bond using a SMARTS pattern.
  2. If not found, then checking for halogen atoms that are attached to another (non‐carbon) atom 
     that in turn is bound to a carbon atom. This “indirect” connection is used to “rescue” cases 
     such as N‐bromosuccinimide or 1‐bromoindole.
Note: In some borderline cases (including some very large natural products) additional rules might be needed.
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    
    A compound is classified as an organohalogen if it contains at least one bond that connects 
    organic matter (a carbon, directly or through a heteroatom linkage) to a halogen atom 
    (F, Cl, Br, or I).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS to catch a direct carbon-halogen bond.
    direct_pattern = Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]")
    if mol.HasSubstructMatch(direct_pattern):
        return True, "Molecule contains a direct carbon–halogen bond"

    # Otherwise, try to “rescue” cases where a halogen is attached to a heteroatom,
    # which in turn is bonded to a carbon. For example, this will catch N‐bromosuccinimide.
    halogen_nums = {9, 17, 35, 53}  # Atomic numbers for F, Cl, Br, I.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in halogen_nums:
            # Look at each neighbor of the halogen.
            for neighbor in atom.GetNeighbors():
                # If neighbor is not carbon, try to see if that neighbor is attached to a carbon.
                if neighbor.GetAtomicNum() != 6:
                    for second_neighbor in neighbor.GetNeighbors():
                        # Ignore the halogen itself.
                        if second_neighbor.GetIdx() == atom.GetIdx():
                            continue
                        if second_neighbor.GetAtomicNum() == 6:
                            return True, ("Molecule contains an indirect carbon–halogen connectivity: "
                                          "halogen attached to a heteroatom that is bonded to carbon")
    # If no matching connectivity was found:
    return False, "No (direct or indirect) carbon–halogen connectivity found"

# Example usage (remove or comment-out these lines if integrating into a larger package):
if __name__ == "__main__":
    # Some test SMILES for demonstration:
    test_smiles = [
         "C[C@H]1C[C@@H](C)Br",  # bromine directly on aliphatic carbon.
         "BrN1C(=O)CCC1=O",      # N-bromosuccinimide (Br attached to N, but N is in a succinimide ring containing C).
         "Cn1ccc2ccccc2c1",      # No halogen.
    ]
    for smi in test_smiles:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")