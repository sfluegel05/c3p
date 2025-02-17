"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python
"""
Classifies: CHEBI:7-hydroxyisoflavones
Definition: A hydroxyisoflavone compound having a hydroxy group at the 7‑position.
This implementation attempts to detect the isoflavone core (benzopyran-4-one or chromen-4-one)
using several SMARTS patterns. Then it searches for a free –OH group attached to an aromatic carbon
in the A‑ring of the core (which, by convention, is the benzene ring fused to the heterocycle).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7‑hydroxyisoflavone based on its SMILES string.
    The process is:
    1. Detect an isoflavone core (the benzopyran-4-one/chromen-4-one skeleton)
       using several alternative SMARTS patterns.
    2. For the first core found, assign the A‑ring atoms (the aromatic benzene ring fused to the pyranone)
       and check for a free hydroxyl (-OH) group directly attached to one of these atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 7‑hydroxyisoflavone, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the input SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to make –OH groups explicit.
    mol = Chem.AddHs(mol)
    
    # Define several SMARTS patterns for the isoflavone (benzopyran-4-one) core.
    # Pattern 1 and 2 represent a clean fused system; Pattern 3 is a relaxed variant
    # that allows extra substituents on the A‑ring.
    core_smarts_patterns = [
        "c1ccc2c(c1)oc(=O)c2",   # pattern 1: typical chromen-4-one skeleton
        "c1ccc2c(c1)occ2=O",      # pattern 2: alternate bonding order for the oxygen
        "c1cc(c(c)c1)-c2oc(=O)cc2"  # pattern 3: expects an explicit connection between an aromatic A‑ring and heterocycle,
                                   # allowing extra substituents on the A‑ring.
    ]
    
    core_found = False
    core_atoms = None
    A_ring_indices = None
    # Loop over the patterns until one matches
    for i, smart in enumerate(core_smarts_patterns):
        core_query = Chem.MolFromSmarts(smart)
        if core_query is None:
            continue  # skip if SMARTS is mis‐defined
        matches = mol.GetSubstructMatches(core_query)
        if matches:
            # Take the first matching substructure for simplicity
            core_atoms = set(matches[0])
            core_found = True
            # Determine which atom indices correspond to the A‑ring.
            # For pattern1 and pattern2 ("c1ccc2c(c1)oc(=O)c2" and "c1ccc2c(c1)occ2=O"),
            # we assume that the first five atoms in the match (indices 0-4) are from the benzene (A‑) ring.
            # For pattern3 ("c1cc(c(c)c1)-c2oc(=O)cc2"), we assume the first six atoms (0-5) form the A‑ring.
            if i < 2:
                A_ring_indices = set(matches[0][:5])
            else:
                A_ring_indices = set(matches[0][:6])
            break  # exit after first successful match
    
    if not core_found:
        return False, "Isoflavone core (benzopyran-4-one skeleton) not detected"
    
    # Look for free hydroxyl groups attached to aromatic carbons.
    # This SMARTS finds an aromatic carbon (in a ring) bonded to an oxygen that has a hydrogen.
    freeOH_query = Chem.MolFromSmarts("[c;R]-[OX2H]")
    if freeOH_query is None:
        return False, "Error in SMARTS definition for free OH group"
    
    freeOH_matches = mol.GetSubstructMatches(freeOH_query)
    if not freeOH_matches:
        return False, "No free hydroxyl (-OH) group attached to an aromatic carbon found"
    
    # Check if any free hydroxyl group is attached to one of the A‑ring atoms of the core
    for match in freeOH_matches:
        aromatic_atom_idx = match[0]  # the aromatic carbon bonded to –OH
        if aromatic_atom_idx in A_ring_indices:
            return True, "Isoflavone core detected with a free hydroxy group on the A-ring (likely 7‑OH)"
    
    return False, "Isoflavone core detected but no free hydroxy group found on the A‑ring (7‑OH missing)"

# Example usage (for testing purposes):
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