"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Any carboxylic acid containing two carboxy groups (dicarboxylic acid)

Our approach:
  1. We parse the molecule and work only on its largest fragment (to remove metal counterions).
  2. We reject molecules that contain any amide bonds (C(=O)N) because many peptides and related molecules
     are false positives for “dicarboxylic acid.”
  3. We also reject molecules that contain common metals (such as Na, K, Ca, Fe, etc.).
  4. We count the number of free carboxy groups by matching SMARTS patterns that capture a carboxyl group –
     either protonated (–COOH) or deprotonated (–COO–). Note that the SMARTS is centered on the carbon in 
     the –C(=O)O fragment.
  5. Only if exactly two unique acid carbons are found do we classify the molecule as a dicarboxylic acid.

Note:
  There is an unavoidable balance here. Some compounds (for example, simple amino acids that contain two free 
  acid groups) may be flagged by a naïve SMARTS match, yet chemists sometimes reserve the term “dicarboxylic acid” 
  for non‐peptidic compounds. In our solution we therefore also filter by (a) rejecting compounds with any amide bond 
  (which flags many peptides) and (b) rejecting molecules that contain any counterions/metals.
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid (exactly two free carboxy acid groups)
    based on its SMILES. 

    The algorithm is as follows:
      1. Parse the SMILES; if invalid, return False.
      2. Add explicit hydrogens.
      3. Split the molecule into fragments and choose the largest (main) fragment (to remove salt/counterions).
      4. Reject if any amide bond (C(=O)N) is detected.
      5. Reject if any metal atoms or common counterions are present.
      6. Count carboxy groups using two SMARTS patterns:
           a. For a protonated acid: "[CX3](=O)[OX2H]" 
           b. For a deprotonated acid: "[CX3](=O)[OX1-]"
         (We consider the matching carbon index as the “tag” for a carboxyl group.)
      7. If exactly two unique acid carbons are found, return True. Otherwise, return False with an explanation.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): (True, explanation) if the molecule is classified as a dicarboxylic acid,
                     (False, explanation) otherwise.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (necessary to see O–H groups)
    mol = Chem.AddHs(mol)
    
    # Remove salts/counterions: consider only the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found in molecule"
    main_frag = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Reject molecules containing amide bonds (C(=O)N) since these are common in peptides.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if amide_pattern is not None and main_frag.HasSubstructMatch(amide_pattern):
        return False, "Molecule contains amide bonds, likely a peptide rather than a simple dicarboxylic acid"
    
    # Reject molecules that contain common metals or counterions.
    # Here we reject if any atom symbol is one of the following:
    metals = {"Na", "K", "Ca", "Fe", "Zn", "Cu", "Mg", "Co", "Mn"}
    for atom in main_frag.GetAtoms():
        if atom.GetSymbol() in metals:
            return False, f"Molecule contains metal ({atom.GetSymbol()}); likely a salt form."

    # Define SMARTS patterns for a carboxyl group.
    # Pattern 1: Protonated carboxylic acid (–COOH)
    acid_prot = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    # Pattern 2: Deprotonated carboxylate (–COO-)
    acid_deprot = Chem.MolFromSmarts("[CX3](=O)[OX1-]")
    
    # We'll use a set to record unique carboxyl carbon atom indices (the matching portion is the carbon).
    acid_carbons = set()
    
    # Find matches for protonated acid groups.
    if acid_prot is not None:
        matches = main_frag.GetSubstructMatches(acid_prot)
        for match in matches:
            # In our SMARTS, match[0] corresponds to the acid carbon.
            acid_carbons.add(match[0])
    
    # Find matches for deprotonated acid groups.
    if acid_deprot is not None:
        matches = main_frag.GetSubstructMatches(acid_deprot)
        for match in matches:
            acid_carbons.add(match[0])
    
    num_acid_groups = len(acid_carbons)
    
    # Decide based on the count.
    if num_acid_groups == 2:
        return True, "Molecule contains exactly two carboxylic acid groups in its main fragment"
    elif num_acid_groups < 2:
        return False, f"Found only {num_acid_groups} carboxylic acid group(s); need exactly two."
    else:
        return False, f"Found {num_acid_groups} carboxylic acid groups; dicarboxylic acid must have exactly two."

# For manual testing under __main__
if __name__ == "__main__":
    examples = {
        "fumarylacetic acid": "OC(=O)CC(=O)\\C=C\\C(O)=O",
        "nocardicin A": "N[C@H](CCOc1ccc(cc1)C(=N\\O)\\C(=O)N[C@H]1CN([C@@H](C(O)=O)c2ccc(O)cc2)C1=O)C(O)=O",
        "meso-tartaric acid": "O[C@@H]([C@@H](O)C(O)=O)C(O)=O",
        "3-ethylmalic acid": "CCC(C(O)C(O)=O)C(O)=O",
        "(4Z)-2-oxohept-4-enedioic acid": "C(=C\\CC(C(O)=O)=O)\\CC(O)=O",
        "2-aminoadipic acid": "NC(CCCC(O)=O)C(O)=O",
        "citraconate(2-)": "C\\C(=C/C([O-])=O)C([O-])=O",  # reported false positive previously
        "potassium L-tartrate": "[K+].[K+].O[C@H]([C@@H](O)C([O-])=O)C([O-])=O",
    }
    for name, smi in examples.items():
        res, reason = is_dicarboxylic_acid(smi)
        print(f"NAME: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")