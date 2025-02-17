"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: Semisynthetic derivative
Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis.
This function uses several heuristics to decide whether a molecule looks like a semisynthetic derivative.
Heuristics include:
  - Must be a valid organic molecule.
  - Should have at least one ring (most natural products have ring systems).
  - Should have at least one chiral center to reflect natural product complexity.
  - Should have a molecular weight in a range common for natural products (here we use 200 - 2000 Da).
  - Optionally, a moderate/high fraction of sp3 hybridized carbons typical of natural products.
Note that these are only approximate rules.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is likely to be a semisynthetic derivative based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets heuristic criteria for a semisynthetic derivative, False otherwise.
        str: Explanation of the reasoning.
    """
    # Parse the SMILES string and get an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is organic (contain at least C and H).
    atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    # Check for carbon (atomic number 6) - natural products are organic.
    if 6 not in atoms:
        return False, "Molecule does not contain carbon atoms and is unlikely organic"

    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 2000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the typical range for semisynthetic derivatives"

    # Check number of rings (most natural products have ring systems).
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 1:
        return False, "No ring system found; semisynthetic derivatives are usually derived from a natural product scaffold with rings"

    # Check for presence of stereocenters (chiral centers are common in natural products).
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 1:
        # This is a clue that the molecule might not be derived from a natural product.
        return False, "No chiral centers found; natural products often contain stereocenters"

    # Check for fraction of sp3 carbons as a surrogate for molecular complexity.
    try:
        frac_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    except Exception:
        frac_sp3 = None

    if frac_sp3 is not None and frac_sp3 < 0.3:
        # A low fraction of sp3 carbons might indicate a simpler, more synthetic-type scaffold.
        return False, f"Fraction of sp3 carbons is low ({frac_sp3:.2f}); may not reflect natural product complexity"
    
    # Optionally, check other factors such as rotatable bonds (if too many, the structure might be too flexible).
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds > 30:
        return False, f"Too many rotatable bonds ({rot_bonds}) for a typical semisynthetic derivative"

    # If passes all tests, we assume the molecule may be a semisynthetic derivative.
    return True, ("Molecule has a ring system, stereocenters, appropriate molecular weight, "
                  "and a reasonable fraction of sp3 carbons consistent with a semisynthetic derivative")

# Example usage:
# result, reason = is_semisynthetic_derivative("CC[C@H]1OC(=O)[C@H](C)...")
# print(result, reason)