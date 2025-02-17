"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide - any lipid carrying one or more hydroperoxy substituents.
The heuristic checks for:
  1. The presence of a hydroperoxy (–OOH) group.
  2. The presence of a long contiguous aliphatic chain (at least 8 non-aromatic, non-ring carbons), which may include unsaturations.
  3. A minimal molecular weight cutoff.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide here is defined heuristically as a molecule that both contains a hydroperoxy (–OOH)
    substituent and a long, unbranched aliphatic chain (at least eight contiguous, non-aromatic and non-cyclic carbons).
    A modest molecular weight cutoff is also used.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a lipid hydroperoxide; False otherwise.
        str: Reason explaining the classification outcome.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the hydroperoxy group.
    # The hydroperoxy group (-OOH) is represented by an -OH bonded to an -O atom.
    # SMARTS: [OX2H]-[OX2] matches an oxygen with an attached hydrogen linked to another oxygen (non-ionized).
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2H]-[OX2]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy (–OOH) group found"
    
    # 2. Check for a long aliphatic chain.
    # Instead of using "CCCCCCCC", we now look for 8 contiguous, non-aromatic carbon atoms.
    # The pattern "[C;!a]" specifies a carbon that is not aromatic.
    # The "~" bond operator allows any bond order (single, double, etc.), capturing unsaturated chains as well.
    # This pattern requires eight such atoms connected in a row.
    chain_pattern = Chem.MolFromSmarts("[C;!a]~[C;!a]~[C;!a]~[C;!a]~[C;!a]~[C;!a]~[C;!a]~[C;!a]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    
    # Look for at least one match where none of the atoms are in a ring.
    found_chain = False
    for match in chain_matches:
        if all(not mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            found_chain = True
            break

    if not found_chain:
        return False, "No long contiguous non-aromatic chain (>=8 carbons) found; molecule may not be lipid-like"
    
    # 3. Check the overall molecular weight.
    # Lipids typically have molecular weight above a certain threshold.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:  # threshold adjusted heuristically to reduce false positives
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical lipid"
    
    return True, "Molecule is classified as a lipid hydroperoxide: contains a hydroperoxy group and a long aliphatic chain"

# Example usage:
if __name__ == "__main__":
    # Example SMILES string (one of those known to be a lipid hydroperoxide)
    test_smiles = "O(O)[C@@H](CCCC(O)=O)/C=C/C=C\\C/C=C\\C/C=C\\CC"  # one example: 5S-HpEPE
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)