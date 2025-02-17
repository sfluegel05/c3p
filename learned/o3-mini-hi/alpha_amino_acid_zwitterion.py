"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
#!/usr/bin/env python
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid–zwitterion obtained by transfer of a proton from the carboxyl 
group to the amino group of any alpha–amino acid; major species at pH 7.3.
This version uses a SMARTS-based recognition of the zwitterionic backbone.
"""

from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines whether a molecule represents an alpha-amino-acid zwitterion.
    The heuristic is to search for a skeletal unit containing:
        a positively charged nitrogen (i.e. [N+;H1,H2,H3])
        directly bound to a carbon (the alpha carbon) that has one or two hydrogens,
        and which is further directly attached to a deprotonated carboxylate group (C(=O)[O-]).
    This pattern corresponds to the canonical zwitterionic backbone for alpha–amino acids.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if an alpha-amino-acid zwitterion pattern is found, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We use explicit hydrogens to ensure the hydrogen counts are accurate.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for the zwitterionic backbone:
    # Pattern explanation:
    #   [N+;H1,H2,H3]  matches any nitrogen with a positive formal charge that has at least one hydrogen.
    #   -               requires a bond between that nitrogen and a carbon.
    #   [C;H1,H2]      matches a carbon atom with one or two hydrogen atoms (as expected in an alpha carbon).
    #   -              requires the carbon to be connected to the carboxylate.
    #   (C(=O)[O-])   matches a carboxylate group (deprotonated C(=O)[O-]).
    zwit_pattern = Chem.MolFromSmarts("[N+;H1,H2,H3]-[C;H1,H2]-(C(=O)[O-])")
    
    if zwit_pattern is None:
        return False, "Error in SMARTS pattern definition"
        
    matches = mol.GetSubstructMatches(zwit_pattern)
    if matches:
        # If one or more matches exist then we classify the molecule as an alpha-amino-acid zwitterion.
        # We report the match (first match) as the basis for classification.
        match = matches[0]
        return True, ("Alpha-amino acid zwitterion pattern found: Nitrogen (atom idx {})—"
                      "alpha carbon (atom idx {})—carboxylate (atom idx corresponding to the carboxylate group)."
                      .format(match[0], match[1]))
    else:
        return False, "No zwitterionic alpha-amino acid backbone (protonated nitrogen linked to an alpha-carbon attached to a deprotonated carboxylate) was found"

# Example usage:
if __name__ == "__main__":
    test_examples = [
        "[C@H]1(CCC[NH+]1C)C([O-])=O",               # N-methylproline zwitterion
        "[C@@H]([C@@H](C([O-])=O)[NH3+])(O)C#C",       # L-beta-ethynylserine zwitterion
        "[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)[O-])=O)NC([C@@H](C3=CC=C(C=C3)OCC[C@H](C([O-])=O)[NH3+])[NH3+])=O",  # nocardicin C dizwitterion
        "CC[C@H](C)[C@H]([NH3+])C([O-])=O",            # L-isoleucine zwitterion
        "CCCC[C@H]([NH3+])C([O-])=O",                  # L-2-aminohexanoic acid zwitterion
        "[NH3+][C@@H](CO)C([O-])=O",                   # L-serine zwitterion
    ]
    
    for smi in test_examples:
        res, msg = is_alpha_amino_acid_zwitterion(smi)
        print("SMILES:", smi)
        print("Result:", res)
        print("Reason:", msg)
        print("---------------------------")