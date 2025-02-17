"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
#!/usr/bin/env python
"""
Classifies: 7-hydroxyisoflavones
Definition: A hydroxyisoflavone compound having a hydroxy group at the 7-position.
This function checks for:
  1. A benzopyran-4-one (isoflavone) core.
  2. A free hydroxyl group appended to one of the aromatic carbons of the A-ring.
The SMARTS used here is "c1cc(O)cc2oc(=O)c(c2)c1" which allows for extra substituents
elsewhere in the molecule while requiring that one aromatic carbon in the fused benzene ring
carries a free –OH group (candidate for the 7-OH).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    It does so by checking for a benzopyran-4-one (isoflavone) core with an attached free
    hydroxyl group on the A ring (as required for the 7-position).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a 7-hydroxyisoflavone, False otherwise.
        str: Reason message for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # (Optional) add hydrogens to help the matching if needed.
    mol_with_H = Chem.AddHs(mol)
    
    # Define a SMARTS for the isoflavone core with an appended free OH on the A ring.
    # The pattern "c1cc(O)cc2oc(=O)c(c2)c1" represents a benzopyran-4-one scaffold
    # in which one of the aromatic carbons of the A ring (the fused benzene) carries a –OH.
    core_smarts = "c1cc(O)cc2oc(=O)c(c2)c1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in SMARTS definition for isoflavone core with OH"
    
    # If the substructure is found, we count it as a match.
    if mol_with_H.HasSubstructMatch(core_query):
        return True, "Isoflavone core with free hydroxy group on the A-ring (7-OH) detected"
    else:
        return False, "Isoflavone core with the required free hydroxy group (7-OH) not found"

# Example usage: (these lines can be commented out for production)
if __name__ == "__main__":
    test_smiles = [
        # 7-hydroxyisoflavones examples:
        "CC(C)=CCc1c(O)cc2occ(-c3ccc(O)cc3O)c(=O)c2c1O",             # luteone
        "CC(C)=CCc1c(O)cc2occ(-c3ccc(O)cc3)c(=O)c2c1O",              # wighteone
        "Oc1ccc2c(c1)occ(-c1ccccc1)c2=O",                            # 7-hydroxyisoflavone
        "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)cc(O)c2c1=O",                # isowighteone
        "COc1c(O)cc(cc1CC=C(C)C)-c1coc2cc(O)ccc2c1=O",                # erylatissin A
        "CC(C)=CCc1c(O)c(O)ccc1-c1coc2cc(O)cc(O)c2c1=O",              # 5,7,3',4'-tetrahydroxy-2'-(3,3-dimethylallyl)isoflavone
        "COc1ccc(cc1O)-c1coc2cc(O)ccc2c1=O",                         # calycosin
        "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",                          # formononetin (should be False: no free A-ring OH)
        "COc1cc2c(cc1O)occ(-c1ccc(O)cc1)c2=O",                       # glycitein
        "OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)cc(O)c2c1occ(-c1ccc(O)cc1)c2=O", # genistein 8-C-glucoside (likely False due to glycosylation)
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
        "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",                            # daidzein (should be False: missing 7-OH in the A-ring)
    ]
    for s in test_smiles:
        result, reason = is_7_hydroxyisoflavones(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")