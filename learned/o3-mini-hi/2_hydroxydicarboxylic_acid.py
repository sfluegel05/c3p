"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid

Definition:
  Any dicarboxylic acid carrying a hydroxy group on the carbon atom
  at position alpha to a carboxyl group.
  
This improved algorithm does the following:
  1. Parses the SMILES (adding explicit hydrogens).
  2. Rejects molecules with nitrogen atoms.
  3. Finds carboxyl (–COOH) groups using the SMARTS "[CX3](=O)[OX2H]" and requires exactly 2.
  4. For each carboxyl group, examines the carbon bonded to the carboxyl carbon (the “alpha-carbon”)
     provided that this candidate is not in a ring or aromatic.
  5. Checks whether that candidate carbon carries an –OH (with at least one hydrogen attached)
     that is not part of the carboxyl group.
     
Note:
  This is a heuristic approach. Some borderline cases (e.g. unsaturated alpha-carbons)
  may be handled differently. The extra requirement for the candidate alpha-carbon 
  not to be in a ring was introduced to help weed out many false positives.
"""

from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    Improved over the previous algorithm by relaxing the sp3 requirement on the alpha carbon
    yet filtering out many cyclic or decorated molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise.
        str : Reason for classification.
    """
    # Parse the SMILES and add hydrogens to help with substructure matching.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Reject any molecule with nitrogen atoms (avoid peptides or N-decorated structures).
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Molecule contains nitrogen atoms which are not expected"
    
    # Define a SMARTS pattern for carboxyl groups: –COOH.
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    if carboxyl_pattern is None:
        return False, "Error creating carboxyl SMARTS pattern"
    
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl groups found"
        
    # Deduplicate based on the carboxyl carbon (the first atom in the match).
    carboxyl_carbons = set(match[0] for match in carboxyl_matches)
    if len(carboxyl_carbons) != 2:
        return False, f"Molecule has {len(carboxyl_carbons)} carboxyl groups (exactly 2 required)"
    
    # Helper: check if an oxygen atom is a genuine hydroxyl (has at least one hydrogen).
    def is_hydroxyl(oxygen_atom):
        if oxygen_atom.GetAtomicNum() != 8:
            return False
        return any(n.GetAtomicNum() == 1 for n in oxygen_atom.GetNeighbors())
    
    # For each carboxyl group (acid carbon), look at its neighbors.
    # The idea: the alpha carbon is the one bonded directly to the carboxyl carbon.
    # To reduce wrong assignments we require that candidate alpha carbons are not in rings and are non‐aromatic.
    for acid_idx in carboxyl_carbons:
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        for neighbor in acid_atom.GetNeighbors():
            # We only consider carbon atoms that are not themselves part of a carboxyl group.
            if neighbor.GetAtomicNum() != 6:
                continue
            if neighbor.GetIdx() in carboxyl_carbons:
                continue
            # Require the candidate alpha carbon to be non–aromatic and not in a ring.
            if neighbor.GetIsAromatic() or neighbor.IsInRing():
                continue
            # Now look for an –OH substituent on the candidate carbon.
            for subnbr in neighbor.GetNeighbors():
                # Skip the bond going back to the acid carbon.
                if subnbr.GetIdx() == acid_idx:
                    continue
                # Check if this neighbor is oxygen and is a hydroxyl (has at least one hydrogen)
                if subnbr.GetAtomicNum() == 8 and is_hydroxyl(subnbr):
                    return True, "Molecule is a 2-hydroxydicarboxylic acid with an alpha hydroxy substituent"
    
    return False, "No suitable alpha-carbon bearing a hydroxyl group adjacent to a carboxyl group was found"


# Uncomment the code below to run tests
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("C[C@H](C(O)=O)[C@@](C)(O)C(O)=O", "(2R,3S)-2,3-dimethylmalic acid"),  # expected True
        ("O[C@@H](CCC(O)=O)C(O)=O", "(S)-2-hydroxyglutaric acid"),  # expected True
        ("O[C@H](CC(O)=O)C(O)=O", "(R)-malic acid"),  # expected True
        ("CCC(C(O)C(O)=O)C(O)=O", "3-ethylmalic acid"),  # expected True
        ("CC(C)[C@@H]([C@@H](O)C(O)=O)C(O)=O", "(2R,3S)-3-isopropylmalic acid"),  # expected True
        ("CC(C)([C@@H](O)C(O)=O)C(O)=O", "(R)-3,3-dimethylmalic acid"),  # expected True
        ("C[C@@](O)(CC(O)=O)C(O)=O", "D-citramalic acid"),  # expected True
        ("OC(CCC(O)=O)C(O)=O", "2-hydroxyglutaric acid"),  # expected True
        ("OC(C(C1CC1=C)C(O)=O)C(O)=O", "2-hydroxy-3-(2-methylidenecyclopropyl)butanedioic acid"),  # expected True
        ("CC[C@@](O)(CC(O)=O)C(O)=O", "(R)-2-ethylmalic acid"),  # expected True
        ("OC(CCCCCC(O)=O)C(O)=O", "2-hydroxyoctanedioic acid"),  # expected True
        ("C(C(CP(=O)(O)[H])(C(O)=O)O)C(O)=O", "2-phosphinomethylmalic acid"),  # expected True
        ("CC(C(O)=O)C(C)(O)C(O)=O", "2,3-dimethylmalic acid"),  # expected True
        ("CC(C)C(C(O)C(O)=O)C(O)=O", "3-isopropylmalic acid"),  # expected True
        ("C[C@](O)(CC(O)=O)C(O)=O", "L-citramalic acid"),  # expected True
        ("O[C@@H](CC(=O)C(O)=O)[C@@H](O)C(O)=O", "5-dehydro-4-deoxy-D-glucaric acid"),  # expected True
        ("C(C(CP(O)=O)C(O)=O)(C(O)=O)O", "phosphinomethylisomalic acid"),  # expected True
        ("OC(=O)\\C=C/C=C(/O)C(O)=O", "(2E,4Z)-2-hydroxymuconic acid"),  # expected True even if unsaturated
        ("OC(=O)\\C=C\\C=C(/O)C(O)=O", "(2Z,4E)-2-hydroxymuconic acid"),  # expected True 
        ("CC(O)(CC(O)=O)C(O)=O", "citramalic acid"),  # expected True
        ("C(C(C(O)=O)O)(CCSC)C(=O)O", "3-(2-methylthioethyl)malic acid"),  # expected True
        ("OC(CCCC(O)=O)C(O)=O", "2-hydroxyadipic acid"),  # expected True
        
        # Some false positives (examples that should NOT be classified)
        ("O[C@H]1C=C(C=C[C@]1(O)C(O)=O)C(O)=O", "(3S,4R)-3,4-dihydroxycyclohexa-1,5-diene-1,4-dicarboxylic acid"),
        ("OC(CC1=CC=C(O)C=C1)(C(O)C(O)=O)C(O)=O", "piscidic acid"),
        ("CSCCCCCC(C(O)C(O)=O)C(O)=O", "3-(5'-Methylthio)pentylmalic acid"),
        ("O[C@H](CC(=O)C(O)=O)C(O)=O", "D-4-hydroxy-2-oxoglutaric acid"),
    ]
    
    for smi, name in test_examples:
        res, reason = is_2_hydroxydicarboxylic_acid(smi)
        print(f"{name}: {res}, Reason: {reason}")