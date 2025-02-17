"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene
Definition: An acyclic (branched or unbranched) hydrocarbon that contains exactly one carbon–carbon double bond
and has the general formula CnH2n.

Improvement: Instead of relying on a regex applied to the string formula, we explicitly count the number of carbon and hydrogen atoms
(using both explicit and implicit hydrogens) after sanitization so that the computed formula is reliable even when stereochemistry or explicit hydrogens are provided.
Also, we verify that the only atoms present are carbon and hydrogen.
Any molecule that fails one of the following criteria is rejected:
  1. Must be a valid molecule.
  2. Must contain only carbon and hydrogen atoms.
  3. Must be acyclic (no rings).
  4. Must contain exactly one carbon–carbon double bond.
  5. Its total atom counts must exactly match the formula CnH2n.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene (acyclic mono‐unsaturated hydrocarbon)
    based on its SMILES string. Specifically, the molecule must:
      1. Be a valid molecule.
      2. Contain only carbon and hydrogen atoms.
      3. Be acyclic.
      4. Have exactly one C=C double bond (between two carbon atoms).
      5. Have a molecular formula of the form CnH2n.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) where the boolean is True if the molecule qualifies as an alkene, and
               the string provides a reason for the classification outcome.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule to correct valences and check consistency.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Remove explicit hydrogens to get a uniform representation – but note that implicit H counts remain available.
    mol = Chem.RemoveHs(mol)
    
    # Check that the molecule contains only carbon (6) and hydrogen (1)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, "Molecule contains atoms other than carbon and hydrogen"
    
    # Verify that the molecule is acyclic (i.e. it has no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected an acyclic structure"
    
    # Count the number of C=C (carbon–carbon double bonds)
    cc_double_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_count += 1
    if cc_double_count != 1:
        return False, f"Molecule has {cc_double_count} C=C bond(s); expected exactly 1"
    
    # Instead of relying on rdMolDescriptors.CalcMolFormula and regex,
    # count the atoms directly based on each atom’s total hydrogen count.
    c_count = 0
    h_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            c_count += 1
            # For carbon, add the number of attached hydrogens (implicit + explicit)
            h_count += atom.GetTotalNumHs()
        elif atom.GetAtomicNum() == 1:
            # In our sanitized, H-removed molecule, explicit hydrogens are gone so we count them as part of GetTotalNumHs() on heavy atoms.
            # This branch should normally not be reached.
            h_count += 1

    if c_count == 0:
        return False, "No carbon atoms found in molecule"
    
    # For an acyclic mono‐unsaturated hydrocarbon, the formula must be exactly CnH2n.
    if h_count != 2 * c_count:
        return False, f"Molecular formula does not match CnH2n; found C{c_count}H{h_count}"
    
    # If all tests pass, the molecule is classified as an alkene.
    return True, f"Molecule is a valid acyclic alkene with 1 C=C bond and formula C{c_count}H{h_count}"