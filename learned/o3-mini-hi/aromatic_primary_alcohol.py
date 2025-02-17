"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: Aromatic Primary Alcohol

Definition:
  An aromatic primary alcohol is any primary alcohol in which the -OH group is attached
  to a carbon that is itself bonded to an atom that is part of an aromatic ring.
  
Approach:
  1. Parse the SMILES string into an RDKit molecule and add explicit hydrogens.
  2. Identify all primary alcohol groups by searching for the substructure "[CH2][OX2H]".
  3. For each match, let the first atom (CH2) be checked:
       • Confirm it has exactly two heavy (non-hydrogen) neighbors.
       • One neighbor must be the hydroxyl oxygen (already part of the match).
       • The other neighbor must be aromatic.
  4. If a match meets these criteria, the molecule qualifies.
  5. Otherwise, return False with an appropriate explanation.
"""

from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
  
    An aromatic primary alcohol is one where the hydroxyl (-OH) group is attached
    to a primary carbon (CH2) that is directly bonded to an aromatic ring.
  
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains an aromatic primary alcohol group, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are accurate.
    mol = Chem.AddHs(mol)
    
    # SMARTS pattern for a primary alcohol group: CH2-OH.
    # This finds a CH2 group directly bonded to an oxygen that is in turn bonded to a hydrogen.
    primary_alcohol_smarts = "[CH2][OX2H]"
    patt = Chem.MolFromSmarts(primary_alcohol_smarts)
    if patt is None:
        return False, "Could not create substructure pattern"
    
    # Look for all substructure matches
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return False, "No primary alcohol (CH2-OH) group found"
    
    # Each match returns a tuple like (idx_CH2, idx_OH)
    for match in matches:
        idx_ch2, idx_oh = match
        ch2_atom = mol.GetAtomWithIdx(idx_ch2)
        oh_atom  = mol.GetAtomWithIdx(idx_oh)
        
        # Count heavy (non-hydrogen) neighbors of the CH2 atom.
        heavy_neighbors = [nbr for nbr in ch2_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        
        # For a primary alcohol, the CH2 should be attached to exactly two heavy atoms:
        # one is the O of the OH group and the other should be from an aromatic ring.
        if len(heavy_neighbors) != 2:
            continue  # This CH2 is not primary.
        
        # Identify the neighbor which is not the oxygen.
        non_oh_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetIdx() != oh_atom.GetIdx()]
        if len(non_oh_neighbors) != 1:
            continue  # Something unexpected.
            
        substituent = non_oh_neighbors[0]
        # Check that this substituent is part of an aromatic ring.
        if not substituent.GetIsAromatic():
            continue  # The non-OH neighbor is not aromatic.
        
        # Also, double-check that the CH2 atom has exactly two hydrogens.
        # (This is because a CH2 should have exactly 2 hydrogens if it is not substituted further.)
        if ch2_atom.GetTotalNumHs() != 2:
            continue  # Not a CH2 group.
        
        # We found a valid aromatic primary alcohol group:
        return True, "Molecule contains an aromatic primary alcohol group"
    
    # If we finish looping through all matches without success, then none met the criteria.
    return False, "No aromatic primary alcohol group found"


# Optional testing block (can be removed if using as a module)
if __name__ == "__main__":
    test_examples = {
        "4-acetoxybenzyl alcohol": "CC(=O)OC1=CC=C(CO)C=C1",
        "3-pyridinemethanol": "C1=CC(=CN=C1)CO",
        "triptohypol A": "O=C1C=C2[C@@]3([C@]([C@]4([C@@](CC3)(CC[C@](C4)(C(=O)O)C)C)[H])(CC[C@]2(C=5C1=C(C(OC)=C(O)C5)CO)C)C)C",
        "(7R,8R)-AGI-B4": "COC(=O)[C@H]1[C@H](O)C=Cc2oc3cc(CO)cc(O)c3c(=O)c12",
        "virensol A": "CCCCC\\C=C\\C=C\\C=C\\C1=C(CO)C(O)=CC=C1O",
        "4-amino-5-hydroxymethyl-2-methylpyrimidine": "Cc1ncc(CO)c(N)n1",
        "tremuloidin": "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](OC(=O)C2=CC=CC=C2)[C@@H]1OC=3C(=CC=CC3)CO)CO",
        "dihydropyriculariol": "C1=CC=C(C(=C1/C=C/C=C/[C@H]([C@H](C)O)O)CO)O",
        "1',4-dihydroxymidazolam": "C1(=NC(C=2N(C=3C1=CC(=CC3)Cl)C(=NC2)CO)O)C=4C=CC=CC4F",
        "2-methylbenzyl alcohol": "CC1=C(CO)C=CC=C1",
        "4-aminopyridine-3-methanol": "C=1C=NC=C(C1N)CO",
        "5-hydroxymethyl-2-furoic acid": "OCC=1OC(C(O)=O)=CC1",
        "1-hydroxymidazolam": "C1(=NCC=2N(C=3C1=CC(=CC3)Cl)C(=NC2)CO)C=4C=CC=CC4F",
        "12-hydroxynevirapine": "O=C1NC=2C(N(C3CC3)C=4N=CC=CC41)=NC=CC2CO",
        "michigazone": "COC1=CC2=NC3=C(OC2=C(OC)C1=O)C=CC(CO)=C3",
        "gentisyl alcohol": "C=1(C=C(C(=CC1)O)CO)O",
        "3,5-dimethylbenzyl alcohol": "CC1=CC(CO)=CC(C)=C1",
        "2-hydroxy-4-hydroxymethylbenzylidenepyruvic acid": "[H]C(=CC(=O)C(O)=O)c1ccc(CO)cc1O",
        "GSK1016790A": "C=1C2=C(C=CC1)SC(=C2)C(=O)N[C@@H](CC(C)C)C(=O)N3CCN(CC3)C(=O)[C@H](CO)NS(C=4C=CC(=CC4Cl)Cl)(=O)=O",
        "5-(hydroxymethyl)cytosine": "Nc1nc(=O)[nH]cc1CO",
        "salicin": "OC[C@H]1O[C@@H](Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O",
        "2-pyrimidinemethanol": "OCC1=NC=CC=N1",
        "(S)-oxamniquine": "CC(C)NC[C@@H]1CCc2cc(CO)c(cc2N1)[N+]([O-])=O",
        "1-butylpyrraline": "CCCCn1c(CO)ccc1C=O",
        "2-hydroxymethyl-3-pentylphenol": "CCCCCC1=C(CO)C(O)=CC=C1",
        "1-neopentylpyrraline": "CC(C)(C)Cn1c(CO)ccc1C=O",
        "1-propylpyrraline": "CCCCn1c(CO)ccc1C=O",
        "1-(5-carboxypentyl)pyrraline": "OCc1ccc(C=O)n1CCCCCC(O)=O",
        "lificiguat": "C(N1N=C(C2=C1C=CC=C2)C=3OC(=CC3)CO)C=4C=CC=CC4",
        "4-(hydroxymethyl)-2-propylfuran-3-carboxylic acid": "CCCC1=C(C(O)=O)C(CO)=CO1",
        "virensol B": "CCCCC\\C=C\\C=C\\C=C\\C1=C(CO)C(O)=CC=C1",
        "2,4-dimethylbenzyl alcohol": "CC1=CC=C(CO)C(C)=C1",
        "(R)-oxamniquine": "CC(C)NC[C@H]1CCc2cc(CO)c(cc2N1)[N+]([O-])=O",
        "3,4-dimethylbenzyl alcohol": "CC1=CC=C(CO)C=C1C",
        "Aloe emodin": "OCc1cc(O)c2C(=O)c3c(O)cccc3C(=O)c2c1",
        "4-cyanobenzyl alcohol": "C1(=CC=C(C=C1)C#N)CO",
    }
    
    for name, smi in test_examples.items():
        res, reason = is_aromatic_primary_alcohol(smi)
        print(f"{name}: {res} - {reason}")