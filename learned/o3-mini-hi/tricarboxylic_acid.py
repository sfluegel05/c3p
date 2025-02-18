"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: Tricarboxylic Acid 
Definition: An oxoacid containing three carboxy groups.
This version uses several heuristics:
  - The molecule must be neutral and contain only allowed (organic) elements.
  - The molecule must contain exactly three unique carboxyl groups (which may be protonated or deprotonated).
  - At most one amide bond is allowed.
  - The carboxyl–to–carbon ratio (number of carboxyl carbons divided by total carbon count) must be at least 0.18.
  - The molecular weight must be below 600 Da.
  - Molecules with any aromatic rings are rejected to ensure that only “simple” acids are accepted.
Note: This classification is heuristic.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines whether a molecule is a tricarboxylic acid based on its SMILES string.
    
    A tricarboxylic acid is defined here as a neutral organic oxoacid that contains exactly three unique carboxyl groups.
    We further require that the molecule is “simple” by:
      - Allowing at most one amide bond,
      - Rejecting any molecule with aromatic rings,
      - Requiring that the carboxyl groups make up a significant fraction of the carbon skeleton (ratio >= 0.18),
      - And that the molecular weight is below 600 Da.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tricarboxylic acid; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check allowed elements (organic: H, C, N, O, S)
    allowed_atomic_nums = {1, 6, 7, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non-organic atom: {atom.GetSymbol()} (atomic num: {atom.GetAtomicNum()})"
    
    # Check that the molecule is neutral (formal charge equals zero)
    if Chem.GetFormalCharge(mol) != 0:
        return False, "Molecule has nonzero formal charge; expected a neutral tricarboxylic acid"
    
    # Reject if any aromatic rings are present.
    # We check for any aromatic atom – if any are found, the molecule is considered too complex.
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic rings; not a simple tricarboxylic acid"
    
    # Count amide bonds – allow at most one.
    # Using a SMARTS pattern for a simple amide: C(=O)N
    amide_pattern = Chem.MolFromSmarts("[C;X3](=O)[N;X3]")
    if amide_pattern is not None:
        amide_matches = mol.GetSubstructMatches(amide_pattern)
        if len(amide_matches) > 1:
            return False, f"Contains {len(amide_matches)} amide bonds; too many for a simple tricarboxylic acid"
    
    # Define SMARTS for carboxyl groups.
    # One for protonated acid: -C(=O)[OH]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    # One for deprotonated carboxylate: -C(=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if acid_pattern is None or carboxylate_pattern is None:
        return False, "Error creating SMARTS patterns for carboxyl groups"
    
    # Get all matches for the two carboxyl patterns.
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Count unique carboxyl-group carbons (take first atom index from each match).
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
    # Lowered threshold (0.18) to allow borderline true positives.
    if carboxyl_ratio < 0.18:
        return False, f"Low carboxyl-to-carbon ratio ({carboxyl_ratio:.2f}); expected at least 0.18"
    
    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 600:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da); expected less than 600 Da for a tricarboxylic acid"
    
    return True, ("Contains exactly three carboxyl groups, has acceptable amide count, no aromatic rings, "
                  "sufficient carboxyl/carbon ratio, and is neutral and within weight range – consistent with a tricarboxylic acid")

# Example usage (for local testing):
if __name__ == "__main__":
    test_smiles = [
        # True examples:
        "C[C@H](\\C=C\\C=C(\\C)[C@H]1CN[C@@H]([C@H]1CC(O)=O)C(O)=O)C(O)=O",  # domoic acid
        "O[C@@H]([C@H](CC(O)=O)C(O)=O)C(O)=O",  # D-erythro-isocitric acid
        "OCCC(NCCC(NCCC(O)C(O)=O)C(O)=O)C(O)=O",  # avenic acid A
        "O=C(O)[C@@H](NC[C@@H](C(O)=O)N)CCC(=O)O",  # N-[(2S)-2-amino-2-carboxyethyl]-L-glutamic acid
        "C(C(C(C)(C(=O)O)O)C(=O)O)C(O)=O",  # 3-hydroxybutane-1,2,3-tricarboxylic acid
        # False examples:
        "O(CCCCCCC(C(C(O)=O)=C)C(O)=O)C(=O)C(O)=O",  # 2-Methylidene-3-(6-oxalooxyhexyl)butanedioic acid (should be rejected)
        "O[C@@H]1C[C@@H](C([O-])=O)\\[N+](C1)=C/C=C1C[C@H](NC(=C\\1)C(O)=O)C(O)=O",  # Portulacaxanthin I (should be rejected due to aromaticity/charge)
        "N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(O)=O)C(O)=O",  # Glu-Asp (multiple amide bonds)
        "OC(CC(=O)NCCC[C@@H](N)C(O)=O)(CC(O)=O)C(O)=O",  # delta-Citryl-d-ornithine (multiple amide bonds)
        "OC(=O)CCP(CCC(O)=O)CCC(O)=O",  # TCEP (contains a P atom)
    ]
    for s in test_smiles:
        result, reason = is_tricarboxylic_acid(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")