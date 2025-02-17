"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer is defined as a mixture, which is composed of macromolecules (usually organic)
of different kinds. In many cases a polymer mixture is represented by a dot-separated SMILES string.
This implementation uses two heuristics:
  (1) If more than one chemically distinct (organic) fragment is present, the mixture is classified as a polymer.
  (2) If only one unique organic fragment is found but the overall mass is not dominated by that fragment (i.e. counter‐ions
      make up >5% of the total mass) and the organic fragment occurs only once, then the substance is considered a polymer.
If no clear decision can be made the function may return (None, None).
Note: This heuristic will not perfectly classify every example.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer mixture.
    
    The approach is:
      1. Parse the SMILES into a molecule and split into disconnected fragments.
      2. Mark fragments as organic if they contain at least one carbon.
      3. For organic fragments, compute a canonical SMILES (to check chemical uniqueness)
         and also compute the molecular weight.
      4. Let total_mass be the summed molecular weight of all fragments.
      5. If two or more unique organic fragments are present then (heuristically) the substance
         is classified as a polymer.
      6. Otherwise, if exactly one unique organic fragment is present but there are other (non‐organic)
         fragments present (i.e. a salt form) then we compute the mass fraction of the organic component.
         If the organic component accounts for less than 95% of the overall mass and it occurs only
         once (rather than a mere multiplicity of the same unit) then classify as polymer.
      7. In all other cases, report that the SMILES does not satisfy our polymer criteria.
      
    Args:
        smiles (str): SMILES string representing the compound.
        
    Returns:
        bool: True if classified as a polymer, False otherwise.
        str: A reason for the decision.
    """
    # Try to parse input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get fragments as separate molecules.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if not fragments:
        return False, "No fragments found in the molecule"
    
    total_mass = sum(rdMolDescriptors.CalcExactMolWt(frag) for frag in fragments)
    
    # We will separate fragments into organic vs inorganic.
    # (Here we call a fragment organic if it contains at least one carbon atom.)
    organic_frags = []
    inorganic_frags = []
    for frag in fragments:
        has_c = any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms())
        if has_c:
            organic_frags.append(frag)
        else:
            inorganic_frags.append(frag)
    
    # Build a dictionary of unique organic fragments.
    # Key: canonical SMILES; Value: list of (fragment, mass)
    unique_org = {}
    total_org_mass = 0.0
    for frag in organic_frags:
        can_smiles = Chem.MolToSmiles(frag, canonical=True)
        mw = rdMolDescriptors.CalcExactMolWt(frag)
        total_org_mass += mw
        if can_smiles in unique_org:
            unique_org[can_smiles]['count'] += 1
            unique_org[can_smiles]['mass'] += mw
        else:
            unique_org[can_smiles] = {'count': 1, 'mass': mw}
    
    n_unique_org = len(unique_org)
    
    # Heuristic decision branch
    if n_unique_org >= 2:
        # At least two chemically distinct organic fragments.
        return True, ("Polymer detected: found at least two distinct organic fragments "
                      "in the mixture (each representing a macromolecular component).")
    elif n_unique_org == 1:
        # Only one unique organic fragment.
        # Many salts have one large organic component and small counterions.
        # If the overall mixture has more than one fragment then check the mass fraction.
        if len(fragments) > 1:
            # Get the mass and occurrence of the unique organic fragment.
            key = list(unique_org.keys())[0]
            org_record = unique_org[key]
            # If the same organic unit occurs by simple replication (multiplicity > 1),
            # we consider it as likely a salt (e.g. quinacetol sulfate).
            if org_record['count'] > 1:
                return False, ("Not classified as polymer: only one unique organic fragment (with duplicates), "
                               "suggesting a salt rather than a polymer mixture.")
            # Otherwise (organic unit occurs once) compare its mass to the total.
            fraction = org_record['mass'] / total_mass
            if fraction < 0.95:
                return True, ("Polymer detected: one large organic fragment accompanied by other species "
                            "reducing its overall mass fraction (<95%), suggesting a polymeric salt formulation.")
            else:
                return False, ("Not classified as polymer: one organic component dominates the mass, "
                               "suggesting a pure compound (or salt) rather than a polymer mixture.")
        else:
            return False, "Not classified as polymer: only one fragment present."
    else:
        # No organic fragments present.
        return False, "Not classified as polymer: no organic (carbon‐containing) fragments found."
    
    # If unable to decide, return (None, None)
    return None, None

# Example usage:
if __name__ == "__main__":
    # You can test a few examples.
    test_smiles_list = [
        # Expected to be polymer by our heuristic
        "C1(=C(C=C2C(=C1OCCN(C)C)C3=C(N2)[C@]4(CC)C)Cl).Cl",  # a simplified mimic of a polymer salt
        "C1CNCCN1.OC(CC(CC(=O)O)(O)C(=O)O).C2CNCCN2.OC(CC(CC(=O)O)(O)C(=O)O",  # mimic of piperazine citrate
        "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(O)=O).CCCCCC",  # mimic of dimethylaminoethyl reserpilinate dihydrochloride
        # Not a polymer:
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",  # one long chain (a macromolecule but not a mixture)
        "C1CCCC1",  # a cyclic compound (only one fragment)
    ]
    
    for smi in test_smiles_list:
        result, reason = is_polymer(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("------")