"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python
"""
Classifies: CHEBI:7-hydroxyisoflavones
Definition: A hydroxyisoflavone compound having a hydroxy group at the 7-position.
This implementation first attempts to detect an isoflavone core (benzopyran-4-one skeleton)
using one of two SMARTS patterns (to be more permissive for additional substituents).
Then it searches for a free hydroxy (-OH) group attached directly to an aromatic carbon 
within the benzene (A-) ring portion of the core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7‑hydroxyisoflavone based on its SMILES string.
    It does so by checking for:
      1. An isoflavone core (benzopyran-4-one skeleton) using alternative SMARTS patterns.
      2. A free hydroxy (–OH) substituent on an aromatic carbon that is part of the benzene ring (A-ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 7‑hydroxyisoflavone, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adding explicit hydrogens so that free OH groups are clearly represented.
    mol = Chem.AddHs(mol)
    
    # Step 1: Search for the isoflavone core (benzopyran-4-one).
    # Use two alternative SMARTS patterns to account for additional substituents.
    core_smarts_patterns = [
        "c1ccc2c(c1)occ2=O",   # original pattern
        "c1ccc2c(c1)oc(=O)c2"  # slightly modified pattern
    ]
    
    core_matches = []
    core_query = None
    for smarts in core_smarts_patterns:
        core_query = Chem.MolFromSmarts(smarts)
        if core_query is None:
            continue  # move on if one of the SMARTS is mis‐defined
        matches = mol.GetSubstructMatches(core_query)
        if matches:
            core_matches = matches[0]  # use the first match for simplicity
            break
    if not core_matches:
        return False, "Isoflavone core (benzopyran-4-one skeleton) not detected"
    
    # For diagnostic purposes, note how many atoms were matched in the core.
    core_atoms = set(core_matches)
    
    # Step 2: Look for a free hydroxy group attached directly to an aromatic carbon.
    # We use the SMARTS "[c;R]-[OX2H]" to look for an aromatic atom bonded to an -OH.
    freeOH_query = Chem.MolFromSmarts("[c;R]-[OX2H]")
    if freeOH_query is None:
        return False, "Error in SMARTS definition for free OH group"
    
    freeOH_matches = mol.GetSubstructMatches(freeOH_query)
    if not freeOH_matches:
        return False, "No free hydroxyl (–OH) group attached to an aromatic carbon found"
    
    # Check if any free hydroxyl is attached to an aromatic carbon that is part of the core.
    # In our approach the core match covers both fused rings.
    # We assume that the A-ring (bearing the 7-OH) is within the atoms of the core.
    for match in freeOH_matches:
        aromatic_atom_idx = match[0]  # aromatic carbon attached to -OH
        # If this aromatic carbon is in the detected core then we consider it as evidence
        # of a free OH on the isoflavone core (likely at the 7-position).
        if aromatic_atom_idx in core_atoms:
            return True, "Isoflavone core detected with a free hydroxy group on the aromatic core (likely 7-OH)"
    
    return False, "Isoflavone core present but free hydroxyl group on the aromatic core (7-OH) not found"


# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "CC(C)=CCc1c(O)cc2occ(-c3ccc(O)cc3O)c(=O)c2c1O",   # luteone
        "CC(C)=CCc1c(O)cc2occ(-c3ccc(O)cc3)c(=O)c2c1O",    # wighteone
        "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O",                 # 7-hydroxyisoflavone
        "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)cc(O)c2c1=O",     # isowighteone
        "COc1c(O)cc(cc1CC=C(C)C)-c1coc2cc(O)ccc2c1=O",     # erylatissin A
        "CC(C)=CCc1c(O)c(O)ccc1-c1coc2cc(O)cc(O)c2c1=O",   # 5,7,3',4'-tetrahydroxy-2'-(3,3-dimethylallyl)isoflavone
        "COc1ccc(cc1O)-c1coc2cc(O)ccc2c1=O",              # calycosin
        "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",               # formononetin (should be False: no free A-ring OH)
        "COc1cc2c(cc1O)occ(-c1ccc(O)cc1)c2=O",            # glycitein
        "OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)cc(O)c2c1occ(-c1ccc(O)cc1)c2=O", # genistein 8-C-glucoside
        "CC(C)(O)CCc1cc(ccc1O)-c1coc2cc(O)cc(O)c2c1=O",   # isowigtheone hydrate
        "CC1(C)Oc2c(O)cc(cc2C=C1)-c1coc2cc(O)cc(O)c2c1=O",  # semilicoisoflavone B
        "Oc1cc(O)c2c(c1)occ(-c1ccc(O)c(O)c1)c2=O",         # orobol
        "COc1ccc(c(O)c1)-c1coc2cc(O)ccc2c1=O",              # 2'-hydroxyformononetin (should be False if OH not on A-ring)
        "CC(C)=CCc1c(O)c(CC=C(C)C)c2occ(-c3ccc(O)cc3)c(=O)c2c1O", # 5,7,4'-trihydroxy-6,8-diprenylisoflavone
        "CC(C)=CCc1c(O)c(CC(O)C(C)=C)c(O)c2c1occ(-c1ccc(O)c(O)c1)c2=O", # millewanin G
        "COC1=C(CC=C(C)C)C(=C(O)C=C1O)C1=COC2=CC(O)=CC=C2C1=O", # kwakhurin
        "Oc1ccc(cc1)-c1coc2cc(O)c(O)cc2c1=O",             # 4',6,7-trihydroxyisoflavone
        "CC(C)=CCc1c(O)ccc(c1O)-c1coc2cc(O)cc(O)c2c1=O",   # licoisoflavone A
        "COc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O",              # 7,4'-dihydroxy-3'-methoxyisoflavone
        "CC(C)=CCc1c(O)ccc(C2COc3cc(O)cc(O)c3C2=O)c1O",    # dihydrolicoisoflavone A
        "COc1cc(OC)c(cc1OC)-c1coc2cc(O)c(OC)cc2c1=O",      # dalnigrein
        "Oc1cc(O)c2c(c1)occ(-c1ccc3OCOc3c1)c2=O",          # 5-hydroxypseudobaptigenin
        "CC1(C)CCc2c(O)ccc(C3COc4cc(O)cc(O)c4C3=O)c2O1",   # glyasperin F
        "OC1=CC=C(C=C1)C1=COC2=C(C(O)=CC(O)=C2)C1=O",       # genistein
        "COc1ccc(cc1)-c1coc2cc(O)cc(O)c2c1=O",             # biochanin A
        "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O",         # neobavaisoflavone
        "C1=C(C(=O)C=2C=CC(O)=C(C2O1)CC=C(C)C)C=3C=CC(O)=CC3", # 8-prenyldaidzein
        "COc1c(O)cc2occ(-c3ccc(O)cc3)c(=O)c2c1O",         # tectorigenin
        "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",                 # daidzein (should be False: missing free OH on the A-ring)
    ]
    for sm in test_smiles:
        result, reason = is_7_hydroxyisoflavones(sm)
        print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")