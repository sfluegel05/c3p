"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthoquinone moiety 
            (a fused bicyclic ring system that is “naphthalene‐like” and
             contains two carbonyl (C=O) groups attached to two of the 10 carbons)
            is substituted by at least one hydroxy group directly attached
            to one of the ring carbons.
            
This implementation uses an explicit SMARTS pattern to search for a 1,4–naphthoquinone core,
and then inspects all atoms in the matched core for external –OH groups.

Note that many natural products have additional fused rings, and the pattern
may miss or wrongly count some substituents when the core is “extended.”
This heuristic therefore may need further refinement for high‐accuracy.
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone.
    
    The function:
    1. Parses the SMILES string.
    2. Searches for a naphthoquinone‐like substructure; here we use a SMARTS pattern
       for a 1,4–naphthoquinone core (an aromatic bicyclic system with two C=O groups).
       (This pattern may not catch other isomers.)
    3. For each match (the candidate naphthoquinone core), inspects every atom in the matched
       set for external substituents. In particular, a valid hydroxy substituent is a single bond
       from a core carbon to an oxygen atom that has at least one hydrogen (as determined by GetTotalNumHs()).
    4. If any candidate core has at least one hydroxy substituent, the molecule is classified 
       as a hydroxynaphthoquinone.
       
    Args:
         smiles (str): SMILES string representation of the molecule.
    
    Returns:
         bool: True if successfully classified as hydroxynaphthoquinone, False otherwise.
         str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for a 1,4–naphthoquinone core.
    # This pattern looks for an aromatic bicyclic system “c1ccc2c(c1)C(=O)C=CC2=O”.
    # Adaptations might be needed for edge–cases.
    naphthoquinone_smarts = "c1ccc2c(c1)C(=O)C=CC2=O"
    core_query = Chem.MolFromSmarts(naphthoquinone_smarts)
    if core_query is None:
        return False, "Error in SMARTS definition"
    
    matches = mol.GetSubstructMatches(core_query)
    if not matches:
        return False, "No naphthoquinone core found"

    # For each naphthoquinone match, check for at least one hydroxy (-OH) substituent 
    # attached directly to one of the core atoms.
    for match in matches:
        core_atom_idxs = set(match)
        hydroxy_count = 0
        
        # We iterate over atoms in the core match.
        for idx in core_atom_idxs:
            atom = mol.GetAtomWithIdx(idx)
            # For each neighbor outside of the core, check if it is an oxygen in an -OH group.
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                # Only consider substituents that are not part of the core match.
                if neighbor.GetIdx() in core_atom_idxs:
                    continue
                # Check if the bond is a single bond, the neighbor atom is oxygen,
                # and the oxygen carries at least one hydrogen.
                if (bond.GetBondType() == Chem.BondType.SINGLE and 
                    neighbor.GetAtomicNum() == 8 and 
                    neighbor.GetTotalNumHs() > 0):
                    hydroxy_count += 1
        if hydroxy_count >= 1:
            return True, (f"Found a naphthoquinone core with {hydroxy_count} hydroxy substituent(s) attached")
    
    # If none of the candidate cores are substituted by an -OH group, then reject.
    return False, "Naphthoquinone core(s) found, but none with at least one hydroxy substituent attached"

# (Optional) Testing examples – uncomment to run some tests:
# if __name__ == '__main__':
#     test_examples = [
#         # True positives examples:
#         ("COC1=C(C)C(=O)c2c(O)cc(OC\\C=C(/C)CCC=C(C)C)cc2C1=O", "7-O-geranyl-2-O,3-dimethylflaviolin"),
#         ("CC(=O)OC(CC=C(C)C)C1=CC(=O)c2c(O)ccc(O)c2C1=O", "Acetylshikonin"),
#         ("Oc1ccc(O)c2C(=O)C=CC(=O)c12", "naphthazarin"),
#         ("Oc1cccc2C(=O)C=CC(=O)c12", "juglone"),
#         ("OC1=CC(=O)c2ccccc2C1=O", "lawsone"),
#         ("Cc1cc(O)c2C(=O)C=CC(=O)c2c1", "Ramentaceone"),
#         ("C1=CC=C(C2=C1C(C=C(C2=O)O)=O)O", "2,8-dihydroxy-1,4-naphthoquinone"),
#         # False negative examples (expected to be classified here if the core is found and -OH detected):
#         ("[C@H](C)([C@@H]([C@@H]([C@H](\\C=C\\O[C@]1(OC=2C(C1=O)=C3C(C(C(=C(C3=O)/C=N/N4CCN(CC4)C)[O-])=O)=C(C2C)[O-])C)OC)C)OC(=O)C)[C@H](O)[C@@H]([C@@H](O)[C@@H](C)/C=C/C=C(/C)C(N)=O)C",
#          "rifampicin para-naphthoquinone carboxamide(2-)"),
#     ]
#     for smi, name in test_examples:
#         result, reason = is_hydroxynaphthoquinone(smi)
#         print(f"NAME: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")