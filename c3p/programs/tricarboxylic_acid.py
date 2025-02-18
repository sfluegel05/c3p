"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: Tricarboxylic Acid 
Definition: An oxoacid containing three carboxy groups.
Heuristics used in this version:
  - The molecule must be neutral and contain only allowed (organic) elements.
  - Molecule must contain exactly three unique carboxyl groups (which may be protonated or deprotonated).
  - There is allowance for at most one amide bond.
  - Molecule must NOT display multiple amino‐acid (peptide) cores.
  - The carboxyl–to–carbon ratio (number of carboxyl carbons divided by total carbon count) must be at least 0.18.
  - The molecular weight must be below 600 Da.
  - Molecules with any aromatic rings are rejected to force “simplicity”.
  - Molecules with unusually high flexibility (rotatable bonds > 12) are rejected.
Note: This classification is heuristic.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines whether a molecule is a tricarboxylic acid based on its SMILES string.
    
    A tricarboxylic acid here is defined as a neutral organic oxoacid that contains exactly three unique carboxyl groups.
    In addition, we require that the compound is “simple” in that:
      - It contains at most one amide bond,
      - It does not contain aromatic rings,
      - It does not appear to be peptide‐like (multiple amino acid units),
      - Its carboxyl groups constitute a significant fraction of the carbon skeleton (ratio >= 0.18),
      - It is not too flexible (rotatable bonds <= 12),
      - And its molecular weight is below 600 Da.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule is classified as a simple tricarboxylic acid; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for allowed elements (organic atoms only: H, C, N, O, S).
    allowed_atomic_nums = {1, 6, 7, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non-organic atom: {atom.GetSymbol()} (atomic num: {atom.GetAtomicNum()})"
    
    # Check that the molecule is neutral (formal charge equals zero)
    if Chem.GetFormalCharge(mol) != 0:
        return False, "Molecule has nonzero formal charge; expected a neutral tricarboxylic acid"
    
    # Reject molecules that contain any aromatic rings.
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic rings; not a simple tricarboxylic acid"
    
    # Count amide bonds using a SMARTS pattern for a simple amide (C(=O)N).
    amide_pattern = Chem.MolFromSmarts("[C;X3](=O)[N;X3]")
    if amide_pattern:
        amide_matches = mol.GetSubstructMatches(amide_pattern)
        if len(amide_matches) > 1:
            return False, f"Contains {len(amide_matches)} amide bonds; too many for a simple tricarboxylic acid"
    
    # Check for peptide‐like fragments.
    # Look for the standard amino acid core: a chiral carbon with an amino group adjacent to a carboxyl group.
    aa_pattern = Chem.MolFromSmarts("[C@H](N)C(=O)O")
    if aa_pattern:
        aa_matches = mol.GetSubstructMatches(aa_pattern)
        if len(aa_matches) >= 2:
            return False, "Molecule appears to be peptide‐like (multiple amino acid moieties)"
    
    # Define SMARTS for a carboxyl group.
    # One for protonated form: -C(=O)[OH]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    # One for deprotonated form: -C(=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if acid_pattern is None or carboxylate_pattern is None:
        return False, "Error creating SMARTS patterns for carboxyl groups"
    
    # Get all matches for the carboxyl groups.
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Count unique carboxyl-group carbon atom indices (use the first atom in each match).
    carboxyl_carbon_indices = set()
    for match in acid_matches:
        carboxyl_carbon_indices.add(match[0])
    for match in carboxylate_matches:
        carboxyl_carbon_indices.add(match[0])
    
    n_carboxyl = len(carboxyl_carbon_indices)
    if n_carboxyl != 3:
        return False, f"Found {n_carboxyl} carboxyl group(s); expected exactly 3"
    
    # Count total carbon atoms.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons == 0:
        return False, "No carbon atoms found"
    
    # Compute the carboxyl-to-carbon ratio.
    carboxyl_ratio = n_carboxyl / total_carbons
    if carboxyl_ratio < 0.18:
        return False, f"Low carboxyl-to-carbon ratio ({carboxyl_ratio:.2f}); expected at least 0.18"
    
    # Check molecular weight (should be below 600 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 600:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da); expected below 600 Da"
    
    # Check rotatable bond count to avoid molecules with overly flexible (long) chains.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 12:
        return False, f"Too many rotatable bonds ({n_rotatable}); not a simple tricarboxylic acid"
    
    return True, ("Contains exactly three carboxyl groups, has acceptable amide count, no aromatic rings, "
                  "no multiple amino acid moieties, sufficient carboxyl/carbon ratio, low rotatable bond count, "
                  "and is neutral and within weight range – consistent with a tricarboxylic acid")

# Example usage (for local testing):
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "C[C@H](\\C=C\\C=C(\\C)[C@H]1CN[C@@H]([C@H]1CC(O)=O)C(O)=O)C(O)=O",  # domoic acid
        "O[C@@H]([C@H](CC(O)=O)C(O)=O)C(O)=O",                           # D-erythro-isocitric acid
        "OCCC(NCCC(NCCC(O)C(O)=O)C(O)=O)C(O)=O",                           # avenic acid A
        "O=C(O)[C@@H](NC[C@@H](C(O)=O)N)CCC(=O)O",                         # N-[(2S)-2-amino-2-carboxyethyl]-L-glutamic acid
        "CCCCCCCCCC[C@H](C(O)=O)[C@@](O)(CC(O)=O)C(O)=O",                  # (2S,3S)-2-hydroxytridecane-1,2,3-tricarboxylic acid
        "C(C(C(C)(C(=O)O)O)C(=O)O)C(O)=O",                                 # 3-hydroxybutane-1,2,3-tricarboxylic acid
        "OC(C(CCC(O)=O)C(O)=O)C(O)=O",                                     # homoisocitric acid
        "[H]C(=CC(=O)C(O)=O)C(CC(O)=O)C(O)=O",                             # 5-oxopent-3-ene-1,2,5-tricarboxylic acid
        "OC(=O)\\C=C(/CC(=O)C(O)=O)C(O)=O",                                # (1E)-4-oxobut-1-ene-1,2,4-tricarboxylic acid
        "OC(=O)CCC(C(O)=O)C(=O)C(O)=O",                                    # 2-oxaloglutaric acid
        "CCC(C(=O)O)(C(=O)O)C(=O)O",                                       # 1,1,1-propanetricarboxylic acid
        "O[C@@H](CN1C[C@@H](O)[C@H]1C(O)=O)[C@H](NCC[C@H](O)C(O)=O)C(O)=O",  # 3-hydroxymugineic acid
        "O[C@H]([C@H](CC(O)=O)C(O)=O)C(O)=O",                              # D-threo-isocitric acid
        "OC(=O)C[C@H](\\C=C\\C(=O)C(O)=O)C(O)=O",                          # (2R,3E)-5-oxopent-3-ene-1,2,5-tricarboxylic acid
        "OC(=O)CCC(=O)N[C@@H](CCCC(=O)C(O)=O)C(O)=O",                     # L-2-succinylamino-6-oxoheptanedioic acid
        # False positives (should be rejected):
        "O(CCCCCCC(C(C(O)=O)=C)C(O)=O)C(=O)C(O)=O",                        # 2-Methylidene-3-(6-oxalooxyhexyl)butanedioic acid
        "N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O",                    # Glu-Asp (dipeptide fragment)
        "OC(CC(=O)NCCC[C@@H](N)C(O)=O)(CC(O)=O)C(O)=O",                     # delta-Citryl-d-ornithine (multiple amide bonds)
        "CC(O)=O.N[C@@H](CCCNC(N)=N)C(=O)N[C@@H](CCC(O)=O)C(O)=O",          # L-arginyl-L-glutamic acid monoacetate (dipeptide fragment)
        "OC(C(O)C(O)=O)(CC(O)=O)C(O)=O",                                    # Hydroxycitric acid (if interpreted as having multiple ambiguous groups)
    ]
    for s in test_smiles:
        result, reason = is_tricarboxylic_acid(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")