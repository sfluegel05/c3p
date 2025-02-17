"""
Classifies: CHEBI:72544 flavonoids
"""
#!/usr/bin/env python3
"""
Classifies: Flavonoids (a superclass including flavonoid, isoflavonoid, chalcone, etc.)
Heuristic improvements include a preliminary removal of sugar moieties (which can inflate carbon counts),
followed by Murcko scaffold determination and then substructure matching.

NOTE: This heuristic is necessarily approximate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def remove_sugars(mol):
    """
    A simple heuristic to remove common sugar units.
    This function looks for substructures matching typical sugar motifs and deletes them.
    (This is not guaranteed to remove all glycosidic parts and may fail in edge cases.)
    """
    # Define a couple of sugar SMARTS patterns (e.g. for pyranoses);
    # these patterns are examples and can be extended.
    sugar_smarts_list = [
        "O[C@H]1[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",  # common glucose/pyranose pattern
        "OC1OC(CO)C(O)C(O)C1O",  # alternate depiction
    ]
    new_mol = Chem.Mol(mol)  # work on a copy
    # Loop over each sugar SMARTS and remove all matches (repeating until no match remains).
    for smarts in sugar_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if not pattern:
            continue
        # Remove all occurrences of the sugar substructure (repeating in case multiple copies are present)
        while new_mol.HasSubstructMatch(pattern):
            new_mol = Chem.DeleteSubstructs(new_mol, pattern)
            # It is important to sanitize after deletion.
            try:
                Chem.SanitizeMol(new_mol)
            except Exception:
                pass
    return new_mol

def is_flavonoids(smiles: str):
    """
    Attempts to determine if a molecule is (putatively) a flavonoid based on its SMILES string.
    The heuristic proceeds as follows:
      1. Validate the SMILES.
      2. Remove likely sugar moieties.
      3. Compute the Murcko scaffold (if possible) of the sugar-stripped molecule.
      4. Verify that the “core” has between about 15 and 21 carbon atoms.
      5. Count aromatic rings (using ring info) – expect 2 (chalcones) or 3 (fused cores) aromatic rings.
      6. Search for one of several SMARTS patterns characteristic of flavonoid substructures;
         these include patterns for flavone, isoflavone, flavanone, flavan, and aurone.
      7. Return True if (a) both the skeletal metrics and a core pattern are satisfied;
         otherwise, return False with an explanation.

    Args:
      smiles (str): SMILES string of the molecule

    Returns:
      bool: True if classified in the flavonoid class, False otherwise.
      str: Explanation for the classification.
    """
    # 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Remove likely sugars to get a better picture of the core.
    mol_nosugar = remove_sugars(mol)
    
    # 3. Attempt to compute the Murcko scaffold from the (sugar-reduced) molecule.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol_nosugar)
    except Exception:
        scaffold = None
    
    # If scaffold computed has at least 15 carbons, use it; if not, fallback on mol_nosugar.
    if scaffold:
        core_carbons = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
        if core_carbons >= 15:
            core = scaffold
        else:
            core = mol_nosugar
    else:
        core = mol_nosugar

    # Count carbons in the selected core.
    core_carbons = sum(1 for atom in core.GetAtoms() if atom.GetAtomicNum() == 6)
    if core_carbons < 15:
        return False, f"Core has only {core_carbons} carbon(s); need at least 15 for flavonoid skeleton."
    if core_carbons > 21:
        return False, f"Core has {core_carbons} carbons (expected between 15 and 21 for typical flavonoid skeleton)."

    # 5. Count aromatic rings in the core.
    ring_info = core.GetRingInfo()
    aromatic_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) >= 6 and all(core.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(tuple(sorted(ring)))
    # Remove duplicate rings
    num_aromatic_rings = len(set(aromatic_rings))
    if num_aromatic_rings not in (2, 3):
        return False, f"Core has {num_aromatic_rings} aromatic ring(s); expected 2 (chalcone) or 3 (fused flavonoid core)."

    # 6. Define SMARTS patterns for flavonoid cores.
    fused_patterns = {
        "flavone":    "c1cc2oc(=O)cc(c2c1)",      # typical flavone/flavonol scaffold
        "isoflavone": "c1ccc2c(c1)oc(=O)c(c2)",    # isoflavone-like core
        "flavanone":  "c1ccc2c(c1)C(=O)[C@H](O)cc2",  # flavanone pattern (includes a chiral center)
        "flavan":     "c1ccc2c(c1)CC(O)cc2",       # flavan pattern (without carbonyl)
        "aurone":     "c1ccc2c(c1)OC(=O)C=CC2"     # aurone pattern (benzofuranone core)
    }
    chalcone_pattern = {
        "chalcone": "c1ccc(cc1)C(=O)C=CC2=CC=CC=C2"
    }
    
    # Select target group based on aromatic ring count.
    if num_aromatic_rings == 3:
        target_patterns = fused_patterns
        target_type = "fused-core"
    elif num_aromatic_rings == 2:
        target_patterns = chalcone_pattern
        target_type = "chalcone"
    else:
        return False, "Unexpected aromatic ring count in core."

    # 7. Attempt to match one of the core SMARTS patterns on the sugar-stripped core or full molecule.
    matched_pattern = None
    for name, smarts in target_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if not pattern:
            continue
        if core.HasSubstructMatch(pattern) or mol_nosugar.HasSubstructMatch(pattern) or mol.HasSubstructMatch(pattern):
            matched_pattern = name
            break

    if matched_pattern:
        return True, f"Matches the {matched_pattern} core pattern ({target_type}) with core containing {core_carbons} carbons and {num_aromatic_rings} aromatic ring(s)."
    
    return False, "No recognized flavonoid core pattern detected in the core or full molecule."

# Optional testing block
if __name__ == "__main__":
    # Some examples taken from the description (note that many of these have sugars)
    test_examples = {
        "4',8-Dimethylgossypetin 3-glucoside": "O1C(C(O)C(O)C(O)C1OC2=C(OC3=C(C2=O)C(O)=CC(O)=C3OC)C4=CC(O)=C(OC)C=C4)CO",
        "Erycristin": "O1C2C(C3=C1C(=C(OC)C=C3)CC=C(C)C)COC4=C2C=C(C(O)=C4)CC=C(C)C",
        "Alternanthin": "O1[C@@](CC(O)[C@H](O)[C@H]1C)(C2=C(O)C3=C(OC(=CC3=O)C4=CC(OC)=C(O)C=C4)C=C2O)[H]",
        "calycosin": "COc1ccc(cc1O)-c1coc2cc(O)ccc2c1=O",
        "wogonin": "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1",
    }
    for name, s in test_examples.items():
        result, reason = is_flavonoids(s)
        print(f"{name}: {result} ({reason})")