"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid class
Definition:
  "A diverse group of pharmacologically active secondary metabolite characteristic to Cannabis plant as well as produced naturally in the body by humans and animals. 
   Cannabinoids contain oxygen as a part of the heterocyclic ring or in the form of various functional groups."
   
Heuristic rules used:
  1. The molecule must contain at least one oxygen atom.
  2. If the molecule contains an aromatic ring with a hydroxyl group (a phenol), it is likely a phytocannabinoid.
  3. Alternatively, if the molecule contains an amide or ester functionality (indicating a polar head group) and a long aliphatic chain (as found in endocannabinoids), 
     then classify it as a cannabinoid.
  
If none of these rules are met then classification returns False with an appropriate reason.
"""

from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is in the cannabinoid class based on heuristics.
    
    Approach:
      - Check that the molecule contains oxygen.
      - Look for an aromatic ring with a hydroxyl group (phenol substructure).
      - If not, look for an amide or ester functionality combined with
        a long aliphatic chain (heuristic for endocannabinoids).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as cannabinoid, False otherwise.
        str: Explanation of the classification decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "Molecule does not contain any oxygen atoms, but cannabinoids require oxygen."
    
    # Rule 1: Look for a phenol substructure (an aromatic ring with an -OH group)
    # SMARTS: aromatic benzene ring with at least one hydroxyl group.
    phenol_smarts = "c1ccc(O)cc1"
    phenol_pattern = Chem.MolFromSmarts(phenol_smarts)
    if phenol_pattern and mol.HasSubstructMatch(phenol_pattern):
        return True, "Molecule contains an aromatic ring with a hydroxyl group (phenol), a common feature of phytocannabinoids."
    
    # Rule 2: Look for amide function (C(=O)N) or ester function (C(=O)O)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    has_amide = mol.HasSubstructMatch(amide_pattern) if amide_pattern else False
    has_ester = mol.HasSubstructMatch(ester_pattern) if ester_pattern else False
    
    # Check for long aliphatic chain. This pattern looks for a chain of at least 7 carbon atoms.
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCC")
    has_long_chain = mol.HasSubstructMatch(long_chain_pattern) if long_chain_pattern else False
    
    if (has_amide or has_ester) and has_long_chain:
        return True, "Molecule has a polar head group (amide/ester) and a long aliphatic chain, features common in endocannabinoids."
    
    # As a further (optional) check, we could inspect rings that contain oxygen.
    # For example, many cannabinoids have heterocyclic rings with oxygen.
    # Here we iterate over rings in the molecule and see if any contains an oxygen atom.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # If the ring size is at least 5 and the ring contains oxygen,
        # assume it might be a heterocycle pertinent to cannabinoid structure.
        if len(ring) >= 5:
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                return True, "Molecule contains a heterocyclic ring with oxygen, a possible cannabinoid feature."
    
    # If none of the above structural features are present, then classify as not a cannabinoid.
    return False, "Molecule does not contain key cannabinoid structural features (phenolic aromatic ring or amide/ester with long chain)."