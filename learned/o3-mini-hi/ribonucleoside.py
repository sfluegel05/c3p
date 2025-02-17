"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: any nucleoside where the sugar component is D-ribose (ribonucleoside)
Definition: A ribonucleoside has a nucleobase linked via a glycosidic bond to a D-ribofuranose unit.
This implementation first excludes compounds with phosphorus (to avoid nucleotides), then searches
for common ribose SMARTS patterns, making sure that at least one atom in the matched sugar is bonded
to an aromatic nitrogen (as a proxy for the nucleobase attachment). Finally, if SMARTS do not match,
a heuristic is used to find a five-membered ring (1 oxygen, 4 carbons) that is directly linked to an aromatic N.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (i.e. a nucleoside with a D-ribose sugar)
    based on its SMILES string.

    Approaches:
      1. Reject molecules containing phosphorus (to avoid nucleotides and phosphorylated species).
      2. Attempt to match common SMARTS patterns of ribose providers and verify that at least one atom 
         in the matched substructure is directly bonded to an aromatic nitrogen (nucleobase signature).
      3. If SMARTS matching fails, search the molecule for any five-membered ring with 1 oxygen & 4 carbons,
         then check for a direct bond from one of the ring carbons to a nitrogen (with aromatic character).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a ribonucleoside, False otherwise.
        str: Explanation/reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if molecule contains phosphorus (P, atomic number 15), likely part of a nucleotide.
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely a nucleotide rather than a nucleoside"
    
    # First approach: use common SMARTS for a ribose substructure.
    # Pattern 1: common unmodified ribofuranose unit.
    ribose_smarts1 = "[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
    pattern1 = Chem.MolFromSmarts(ribose_smarts1)
    if pattern1 is None:
        return False, "Error creating SMARTS pattern for ribose (pattern1)"
    # Pattern 2: variant (e.g., for 2'-O-methyl nucleosides).
    ribose_smarts2 = "CO[C@@H]1[C@H](O)[C@@H](CO)O[C@H]1"
    pattern2 = Chem.MolFromSmarts(ribose_smarts2)
    if pattern2 is None:
        return False, "Error creating SMARTS pattern for ribose (pattern2)"
    
    # Helper function to check if a substruct match is linked to an aromatic nitrogen (nucleobase).
    def sugar_attached_to_nucleobase(match):
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider neighbors not in the sugar match.
                if nbr.GetIdx() not in match and nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic():
                    return True
        return False

    # Check pattern 1.
    matches1 = mol.GetSubstructMatches(pattern1)
    for match in matches1:
        if sugar_attached_to_nucleobase(match):
            return True, "Matched common D-ribose pattern (unmodified ribofuranose) attached to a nucleobase"
    # Check pattern 2.
    matches2 = mol.GetSubstructMatches(pattern2)
    for match in matches2:
        if sugar_attached_to_nucleobase(match):
            return True, "Matched common D-ribose pattern with 2'-O-methyl modification attached to a nucleobase"
    
    # Second approach: heuristic search for five-membered ribofuranose rings.
    # Use an explicit hydrogen-added molecule for better perception.
    mol_with_H = Chem.AddHs(mol)
    ring_info = mol_with_H.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Count the number of oxygen and carbon atoms.
            oxygen_idxs = []
            carbon_idxs = []
            for idx in ring:
                atom = mol_with_H.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxygen_idxs.append(idx)
                elif atom.GetAtomicNum() == 6:
                    carbon_idxs.append(idx)
            # Ribofuranose: exactly one oxygen, four carbons.
            if len(oxygen_idxs) == 1 and len(carbon_idxs) == 4:
                # For each carbon in the ring, check if it is bonded to an external aromatic nitrogen.
                for idx in carbon_idxs:
                    atom = mol_with_H.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        # Verify the neighboring atom is outside the ring.
                        if nbr.GetIdx() not in ring and nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic():
                            return True, "Found five-membered ribofuranose ring attached to a nucleobase nitrogen"
    return False, "No ribose moiety attached to a nucleobase found"

# Example usage:
if __name__ == "__main__":
    # Test a few examples:
    examples = {
        "1-methyladenosine": "Cn1cnc2n(cnc2c1=N)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
        "3,4-dihydrozebularine": "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=CCNC1=O",
        "nucleocidin": "NC1=C2N=CN([C@@H]3O[C@](F)(COS(N)(=O)=O)[C@@H](O)[C@H]3O)C2=NC=N1",
        "cytidine": "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1",
        "5'-deoxyadenosine": "C[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",  # should be false as deoxy sugar
        "ATP (should be rejected)": "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O"
    }
    for name, smi in examples.items():
        result, reason = is_ribonucleoside(smi)
        print(f"{name}: {result} ({reason})")