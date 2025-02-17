"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen Compound
A compound containing at least one carbon–halogen bond (where X is F, Cl, Br or I).

This implementation first looks for a direct bond between any carbon and a halogen.
If none is found, it performs an additional search for halogen atoms that are attached 
to any heteroatom (e.g. N, O, S, etc.) that in turn is bound to a carbon.
This is meant to “rescue” cases such as N‐bromosuccinimide or 1‐bromoindole.
Note: In some borderline cases (including some very large natural product molecules)
additional context‐dependent rules could be added.
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    
    A compound is classified as an organohalogen if it contains at least one carbon–halogen bond.
    This function first checks for a direct bond between any carbon (aliphatic or aromatic)
    and a halogen (F, Cl, Br, I). If not found, it then checks for an “indirect” connection,
    where a halogen atom is attached to some heteroatom (non‐carbon, non‐halogen) that itself 
    is bound to a carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. First pass: Look for a direct bond between carbon and halogen.
    # The SMARTS "[#6]-[F,Cl,Br,I]" catches a bond between any carbon (atomic number 6)
    # and any halogen (F, Cl, Br, I). This pattern works for both aromatic and aliphatic atoms.
    direct_pattern = Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]")
    if mol.HasSubstructMatch(direct_pattern):
        return True, "Molecule contains a direct carbon–halogen bond"
    
    # --- 2. Second pass: Look for halogens that are indirectly connected.
    # Sometimes a halogen is bound not directly to carbon but to a heteroatom (e.g. N, O, S),
    # which is itself bound to a carbon. This loop looks for any halogen atom and then 
    # confirms whether any neighbor (which is not carbon) has at least one neighboring carbon.
    halogen_atomic_nums = {9, 17, 35, 53}  # F, Cl, Br, I respectively.
    
    # Iterate over all atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in halogen_atomic_nums:
            # Examine each neighbor of the halogen.
            for neigh in atom.GetNeighbors():
                # If the neighbor is carbon, we would have caught it already.
                # But if not carbon (and not another halogen), then look at its neighbors.
                if neigh.GetAtomicNum() != 6 and neigh.GetAtomicNum() not in halogen_atomic_nums:
                    for neigh2 in neigh.GetNeighbors():
                        if neigh2.GetIdx() == atom.GetIdx():
                            continue  # ignore the original halogen
                        if neigh2.GetAtomicNum() == 6:
                            return True, ("Molecule contains an indirect carbon–halogen connectivity: "
                                          "halogen attached to a heteroatom that is bonded to carbon")
    # If no suitable connectivity is found, then we consider this molecule negative.
    return False, "No (direct or indirect) carbon–halogen connectivity found"


# Example usage: (Comment out or remove when integrating into a larger package)
if __name__ == "__main__":
    # Some test examples including those expected to turn up as organohalogens.
    test_smiles = [
         "C[C@H]1C[C@@H](C)Br",    # direct C–Br bond
         "BrN1C(=O)CCC1=O",         # N‐bromosuccinimide (halogen indirectly attached through N)
         "OC\\C=C\\Cl",            # trans-3-chloroprop-2-en-1-ol, direct C–Cl bond
         "Oc1cccc(F)c1O",          # 3-fluorocatechol, aromatic C–F bond
         "O=c1cc(oc2c(cccc12)-c1ccccc1)N1CCOCC1",  # LY294002; no halogen present
    ]
    for smi in test_smiles:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")