"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python
"""
Classifies: CHEBI:7-hydroxyisoflavones
Definition: A hydroxyisoflavone compound having a hydroxy group at the 7-position.

This function works in two steps:
  1. It searches for an isoflavone (benzopyran-4-one) core.
     We use a looser SMARTS pattern "c1ccc2c(c1)occ2=O" so that extra substituents
     on the ring (such as prenyl groups) do not prevent a match.
  2. It then looks for a free hydroxyl group (–OH) attached directly to one of the
     aromatic carbons of the fused benzene ring in the isoflavone core.
  
If both are found, it is classified as a 7-hydroxyisoflavone.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7‑hydroxyisoflavone based on its SMILES string.
    It does so by checking for an isoflavone core (benzopyran-4-one skeleton)
    plus a free hydroxy (–OH) substituent attached to one of the aromatic carbons of
    the A-ring (which typically corresponds to the 7‑position).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a 7‑hydroxyisoflavone, False otherwise.
        str: Reason message for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adding hydrogens to ensure that our free –OH groups are represented explicitly.
    mol = Chem.AddHs(mol)
    
    # Step 1: Search for the isoflavone core.
    # We use a SMARTS that defines a benzopyran-4-one scaffold.
    # The pattern "c1ccc2c(c1)occ2=O" captures a fused aromatic ring (A-ring) with a heterocycle with a ketone.
    core_smarts = "c1ccc2c(c1)occ2=O"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in SMARTS definition for isoflavone core"
    
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Isoflavone core (benzopyran-4-one skeleton) not detected"
    
    # For simplicity we take the atoms involved in the first found core match.
    # (This match is a tuple of atom indices corresponding to the core SMARTS.)
    core_atoms = set(core_matches[0])
    
    # Step 2: Look for a free hydroxy group attached to an aromatic carbon in the core.
    # Free hydroxyl groups are represented as an oxygen (OX2) having one hydrogen.
    # The SMARTS "[c;R]-[OX2H]" looks for an aromatic (ring) carbon bonded to a free –OH.
    freeOH_smarts = "[c;R]-[OX2H]"
    freeOH_query = Chem.MolFromSmarts(freeOH_smarts)
    if freeOH_query is None:
        return False, "Error in SMARTS definition for free OH group"

    freeOH_matches = mol.GetSubstructMatches(freeOH_query)
    # Now check if any free –OH is attached to an aromatic carbon that is part of the core.
    for match in freeOH_matches:
        # match is a tuple where the first atom is the aromatic carbon neighbor of –OH
        aromatic_c_idx = match[0]
        if aromatic_c_idx in core_atoms:
            return True, "Isoflavone core detected with a free hydroxy group attached to the aromatic core (likely 7-OH)"
    
    return False, "Isoflavone core present but free hydroxy group on the aromatic A-ring (7-OH) not found"

# Example usage (these lines may be removed or commented out in production):
if __name__ == "__main__":
    test_smiles = [
        "CC(C)=CCc1c(O)cc2occ(-c3ccc(O)cc3O)c(=O)c2c1O",             # luteone
        "CC(C)=CCc1c(O)cc2occ(-c3ccc(O)cc3)c(=O)c2c1O",              # wighteone
        "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O",                            # 7-hydroxyisoflavone
        "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)cc(O)c2c1=O",                # isowighteone
        "COc1c(O)cc(cc1CC=C(C)C)-c1coc2cc(O)ccc2c1=O",                # erylatissin A
        "CC(C)=CCc1c(O)c(O)ccc1-c1coc2cc(O)cc(O)c2c1=O",              # 5,7,3',4'-tetrahydroxy-2'-(3,3-dimethylallyl)isoflavone
        "COc1ccc(cc1O)-c1coc2cc(O)ccc2c1=O",                         # calycosin
        "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",                          # formononetin (should be False: no free A-ring OH)
        "COc1cc2c(cc1O)occ(-c1ccc(O)cc1)c2=O",                       # glycitein
        "OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)cc(O)c2c1occ(-c1ccc(O)cc1)c2=O", # genistein 8-C-glucoside
        "CC(C)(O)CCc1cc(ccc1O)-c1coc2cc(O)cc(O)c2c1=O",             # isowigtheone hydrate
        "CC1(C)Oc2c(O)cc(cc2C=C1)-c1coc2cc(O)cc(O)c2c1=O",            # semilicoisoflavone B
        "Oc1cc(O)c2c(c1)occ(-c1ccc(O)c(O)c1)c2=O",                   # orobol
        "COc1ccc(c(O)c1)-c1coc2cc(O)ccc2c1=O",                        # 2'-hydroxyformononetin (should be False if OH is not on the A ring)
        "CC(C)=CCc1c(O)c(CC=C(C)C)c2occ(-c3ccc(O)cc3)c(=O)c2c1O",     # 5,7,4'-trihydroxy-6,8-diprenylisoflavone
        "CC(C)=CCc1c(O)c(CC(O)C(C)=C)c(O)c2c1occ(-c1ccc(O)c(O)c1)c2=O",# millewanin G
        "COC1=C(CC=C(C)C)C(=C(O)C=C1O)C1=COC2=CC(O)=CC=C2C1=O",       # kwakhurin
        "Oc1ccc(cc1)-c1coc2cc(O)c(O)cc2c1=O",                        # 4',6,7-trihydroxyisoflavone
        "CC(C)=CCc1c(O)ccc(c1O)-c1coc2cc(O)cc(O)c2c1=O",             # licoisoflavone A
        "COc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O",                         # 7,4'-dihydroxy-3'-methoxyisoflavone
        "CC(C)=CCc1c(O)ccc(C2COc3cc(O)cc(O)c3C2=O)c1O",              # dihydrolicoisoflavone A
        "COc1cc(OC)c(cc1OC)-c1coc2cc(O)c(OC)cc2c1=O",                # dalnigrein
        "Oc1cc(O)c2c(c1)occ(-c1ccc3OCOc3c1)c2=O",                    # 5-hydroxypseudobaptigenin
        "CC1(C)CCc2c(O)ccc(C3COc4cc(O)cc(O)c4C3=O)c2O1",             # glyasperin F
        "OC1=CC=C(C=C1)C1=COC2=C(C(O)=CC(O)=C2)C1=O",                 # genistein
        "COc1ccc(cc1)-c1coc2cc(O)cc(O)c2c1=O",                        # biochanin A
        "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O",                   # neobavaisoflavone
        "C1=C(C(=O)C=2C=CC(O)=C(C2O1)CC=C(C)C)C=3C=CC(O)=CC3",         # 8-prenyldaidzein
        "COc1c(O)cc2occ(-c3ccc(O)cc3)c(=O)c2c1O",                    # tectorigenin
        "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",                            # daidzein (should be False: missing free OH on the A-ring)
    ]
    for s in test_smiles:
        result, reason = is_7_hydroxyisoflavones(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")