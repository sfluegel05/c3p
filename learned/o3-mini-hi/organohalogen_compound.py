"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen Compound
A compound is considered an organohalogen if it contains at least one carbon–halogen bond
(where the halogen is F, Cl, Br, or I).
This implementation first cleans the molecule by extracting the largest fragment, then
tries to locate a direct carbon–halogen bond using a SMARTS pattern. If none is found,
it attempts to locate an “indirect” connectivity (e.g. a halogen bonded to a heteroatom
that in turn is bonded to a carbon) via two complementary SMARTS patterns.
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    
    A compound is classified as an organohalogen if it contains at least one carbon–halogen bond.
    This function looks both for direct bonds between carbon and a halogen (F, Cl, Br or I),
    and for indirect connectivity where a halogen is attached to a non-carbon (and non-hydrogen)
    atom that is itself bonded to a carbon.
    
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
    
    # To avoid counting small counterion fragments, use only the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    largest = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # --- 1. Direct connectivity: look for a bond between any carbon and a halogen.
    # This SMARTS looks for any carbon (atomic number 6) directly bonded to any of F, Cl, Br, I.
    direct_smarts = "[#6]-[F,Cl,Br,I]"
    direct_pattern = Chem.MolFromSmarts(direct_smarts)
    if largest.HasSubstructMatch(direct_pattern):
        return True, "Molecule contains a direct carbon–halogen bond"
    
    # --- 2. Indirect connectivity: sometimes the halogen is attached to a heteroatom (not C or H)
    # that in turn is bonded to a carbon. Use two complementary SMARTS patterns.
    # Pattern 1: carbon attached to a heteroatom, which is then attached to a halogen.
    indirect_smarts = "[#6]-[!#6&!#1]-[F,Cl,Br,I]"
    indirect_pattern = Chem.MolFromSmarts(indirect_smarts)
    if largest.HasSubstructMatch(indirect_pattern):
        return True, ("Molecule contains an indirect carbon–halogen connectivity via "
                      "a heteroatom bound to carbon")
    
    # Pattern 2: halogen attached to a heteroatom, which is then attached to a carbon.
    indirect_smarts_rev = "[F,Cl,Br,I]-[!#6&!#1]-[#6]"
    indirect_pattern_rev = Chem.MolFromSmarts(indirect_smarts_rev)
    if largest.HasSubstructMatch(indirect_pattern_rev):
        return True, ("Molecule contains an indirect carbon–halogen connectivity via "
                      "a heteroatom bound to carbon")
    
    # If no pattern is matched, then the molecule lacks any carbon–halogen bond.
    return False, "No (direct or indirect) carbon–halogen connectivity found"

# Example usage (for testing; remove or comment out if integrating into a larger package)
if __name__ == "__main__":
    test_smiles = [
         "C[C@H]1C[C@@H](C)Br",     # direct C–Br bond
         "BrN1C(=O)CCC1=O",          # N-bromosuccinimide (halogen indirectly attached via N)
         "OC\\C=C\\Cl",             # trans-3-chloroprop-2-en-1-ol, direct C–Cl bond
         "Oc1cccc(F)c1O",           # 3-fluorocatechol, aromatic C–F bond
         "O=c1cc(oc2c(cccc12)-c1ccccc1)N1CCOCC1",  # LY294002; no halogen connectivity
    ]
    for smi in test_smiles:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")