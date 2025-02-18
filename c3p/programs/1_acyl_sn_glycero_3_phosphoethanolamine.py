"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
This program uses substructure SMARTS to check for the phosphoethanolamine headgroup,
the ester linkage (indicating the acyl chain at the sn-1 position), and finally uses
chiral center assignment to ensure that the molecule has the (R)-configuration.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    
    The molecule must have:
      - A phosphoethanolamine headgroup (as indicated by a phosphate bonded to an ethanolamine fragment).
      - An ester linkage indicating the 1-O-acyl chain attached to the glycerol backbone.
      - A chiral center (on the glycerol moiety) that is assigned the (R) configuration.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified in the target class, False otherwise.
        str: A message explaining the reasoning.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure aromaticity and assign stereochemistry if needed
    Chem.SanitizeMol(mol)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # 1. Check for the phosphoethanolamine headgroup.
    # This SMARTS checks for a phosphorus with one double-bonded oxygen,
    # one hydroxyl oxygen and one oxygen attached to an ethanolamine (OCCN).
    phos_ethanolamine_smarts = "P(=O)(O)OCCN"
    phos_ethanolamine_pattern = Chem.MolFromSmarts(phos_ethanolamine_smarts)
    if not mol.HasSubstructMatch(phos_ethanolamine_pattern):
        return False, "Phosphoethanolamine headgroup (P(=O)(O)OCCN) not found"

    # 2. Check for the ester linkage that attaches a fatty acyl chain to the glycerol.
    # We consider two possibilities for the chiral specification on the glycerol oxygen:
    ester_smarts_1 = "C(=O)O[C@H]"   # using @H for one possible chirality mark
    ester_smarts_2 = "C(=O)O[C@@H]"  # alternate chiral marker
    ester_pattern1 = Chem.MolFromSmarts(ester_smarts_1)
    ester_pattern2 = Chem.MolFromSmarts(ester_smarts_2)
    if not (mol.HasSubstructMatch(ester_pattern1) or mol.HasSubstructMatch(ester_pattern2)):
        return False, "Ester group linking acyl chain to glycerol (C(=O)O[C@H] or C(=O)O[C@@H]) not found"
    
    # 3. Check for a chiral centre with assigned (R) configuration.
    # We use RDKit's chiral center finder.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if not chiral_centers:
        return False, "No defined chiral centers found"
    
    # Look for at least one chiral centre that is assigned as (R).
    r_center_found = False
    for center in chiral_centers:
        idx, config = center  # config is a string "R" or "S"
        if config == "R":
            r_center_found = True
            break
    if not r_center_found:
        return False, "No (R)-configured chiral center found; required (R)-configuration missing"

    # If all checks pass, return True.
    return True, "Molecule contains phosphoethanolamine headgroup, an acyl ester linkage on glycerol, and an (R)-configured chiral center"

# Optional simple self-test block (can be removed if not desired)
if __name__ == "__main__":
    # Example: 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"
    result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(test_smiles)
    print("Test molecule classification:", result)
    print("Reason:", reason)