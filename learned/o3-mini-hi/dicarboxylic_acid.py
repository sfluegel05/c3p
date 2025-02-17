"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Any carboxylic acid containing two carboxy groups (dicarboxylic acid)

Improved Approach:
  1. Parse the SMILES string and add explicit hydrogens.
  2. Consider only the largest fragment (to remove salts/counterions).
  3. Reject molecules that contain common metals.
  4. Identify all carboxyl groups – both protonated ([CX3](=O)[OX2H]) and deprotonated ([CX3](=O)[OX1-]).
     (Here we “tag” a carboxyl group by its carbon atom.)
  5. Only if exactly two unique acid carbons are found do we proceed.
  6. We then compute the shortest path (in bonds) between the two acid carbons. Since in classical
     “dicarboxylic acids” the carboxyl groups are connected by a simple aliphatic (or aromatic) carbon–based 
     backbone, we require that (a) every intermediate atom along the path is a carbon, and (b) the overall path
     length falls within a reasonable range.
  7. If these criteria do not hold then we reject the molecule.
  
Note:
  We no longer outright reject molecules that contain amide bonds because some target compounds (e.g. nocardicin A, 
  methotrexate) do include amide functionality yet are considered dicarboxylic acids by the “exactly two” rule.
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid (exactly two free carboxy acid groups)
    based on its SMILES string.

    The algorithm is as follows:
      1. Parse the SMILES and add explicit hydrogens.
      2. Split the molecule into fragments and choose the largest (main fragment).
      3. Reject if any metal atoms are present.
      4. Identify carboxyl groups using two SMARTS patterns, and record the matching carbon atom.
      5. If exactly two unique acid carbons are found, compute the shortest bond path between them.
         Then require that:
           a. All intermediate atoms on this path are carbon.
           b. The path length (number of bonds) is neither too short (<2) nor too long (>12).
      6. If these checks pass, return True; otherwise, return False with an explanation.

    Args:
        smiles (str): SMILES string representation of the molecule.

    Returns:
        (bool, str): Tuple (True, explanation) for a dicarboxylic acid; else (False, explanation).
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (so that protons on -OH are visible).
    mol = Chem.AddHs(mol)
    
    # Remove salts/counterions: consider only the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found in molecule"
    main_frag = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Reject if any metal atoms or common counterions are present.
    metals = {"Na", "K", "Ca", "Fe", "Zn", "Cu", "Mg", "Co", "Mn"}
    for atom in main_frag.GetAtoms():
        if atom.GetSymbol() in metals:
            return False, f"Molecule contains metal ({atom.GetSymbol()}); likely a salt form."
    
    # Define SMARTS patterns for a carboxyl group.
    acid_prot = Chem.MolFromSmarts("[CX3](=O)[OX2H]")   # protonated acid (–COOH)
    acid_deprot = Chem.MolFromSmarts("[CX3](=O)[OX1-]")   # deprotonated acid (–COO-)
    
    # Use a set to record unique carboxyl carbon atom indices.
    acid_carbons = set()
    
    # Find matches for protonated acid groups.
    if acid_prot is not None:
        for match in main_frag.GetSubstructMatches(acid_prot):
            # In our SMARTS, the matching atom 0 is the carbon of the –COOH group.
            acid_carbons.add(match[0])
    
    # Find matches for deprotonated acid groups.
    if acid_deprot is not None:
        for match in main_frag.GetSubstructMatches(acid_deprot):
            acid_carbons.add(match[0])
    
    num_acid = len(acid_carbons)
    if num_acid != 2:
        if num_acid < 2:
            return False, f"Found only {num_acid} carboxylic acid group(s); need exactly two."
        else:
            return False, f"Found {num_acid} carboxylic acid groups; dicarboxylic acid must have exactly two."
    
    # At this point, exactly two unique acid carbons were found.
    acid_list = list(acid_carbons)
    # Get the shortest path (list of atom indices) between the two acid carbons.
    path = Chem.GetShortestPath(main_frag, acid_list[0], acid_list[1])
    if not path or len(path) < 2:
        return False, "Could not determine a valid connection between the two acid groups."
    
    path_length = len(path) - 1  # number of bonds
    # Impose a reasonable path length
    if path_length < 2:
        return False, "Acid groups appear to be too close (not connected by an intervening chain)."
    if path_length > 12:
        return False, "Acid groups are separated by a long and complex structure; not a simple dicarboxylic acid."
    
    # Check that all intermediate atoms (excluding the acid carbons themselves) are carbon.
    # (In a typical dicarboxylic acid, the two acid groups are connected by a hydrocarbon chain or an aromatic ring.)
    for idx in path[1:-1]:
        atom = main_frag.GetAtomWithIdx(idx)
        if atom.GetSymbol() != "C":
            return False, "Acid groups are not connected by a simple carbon chain (found heteroatoms in the connecting path)."
    
    return True, "Molecule contains exactly two carboxylic acid groups connected by a simple carbon backbone in its main fragment"

# For manual testing
if __name__ == "__main__":
    examples = {
        "fumarylacetic acid": "OC(=O)CC(=O)\\C=C\\C(O)=O",
        "nocardicin A": "N[C@H](CCOc1ccc(cc1)C(=N\\O)\\C(=O)N[C@H]1CN([C@@H](C(O)=O)c2ccc(O)cc2)C1=O)C(O)=O",
        "meso-tartaric acid": "O[C@@H]([C@@H](O)C(O)=O)C(O)=O",
        "3-ethylmalic acid": "CCC(C(O)C(O)=O)C(O)=O",
        "(4Z)-2-oxohept-4-enedioic acid": "C(=C\\CC(C(O)=O)=O)\\CC(O)=O",
        "2-aminoadipic acid": "NC(CCCC(O)=O)C(O)=O",
        "citraconate(2-)": "C\\C(=C/C([O-])=O)C([O-])=O",
        "potassium L-tartrate": "[K+].[K+].O[C@H]([C@@H](O)C([O-])=O)C([O-])=O",
    }
    for name, smi in examples.items():
        res, reason = is_dicarboxylic_acid(smi)
        print(f"NAME: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")