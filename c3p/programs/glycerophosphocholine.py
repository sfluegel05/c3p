"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
Definition: The glycerol phosphate ester of a phosphocholine.
A glycerophosphocholine molecule should include (i) a phosphate group bonded to a glycerol or glycerol-like moiety and
(ii) a choline headgroup (–OCC[N+](C)(C)C) attached via that phosphate.
Note: Many natural examples (e.g., phosphatidylcholines) contain additional acyl chains.
"""

from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine (GPC) based on its SMILES string.
    A glycerophosphocholine is essentially a glycerol phosphate ester of phosphocholine.
    
    Key criteria (not exhaustive):
      1. The molecule must contain a phosphorus atom.
      2. It must have a phosphocholine headgroup, i.e. a phosphate group linked to an oxygen that is in turn bonded
         to an ethylene chain ending in a trimethylammonium (choline) group: OCC[N+](C)(C)C.
         
    Due to variability in how the phosphate is drawn (protonated vs deprotonated), we will try several SMARTS
    patterns for the phosphocholine headgroup.

    Args:
        smiles (str): SMILES string of the molecule to classify

    Returns:
        bool: True if the molecule is classified as glycerophosphocholine, False otherwise.
        str: Explanation of the result.
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, check that the molecule contains at least one phosphorus atom (atomic number 15).
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found; not a phospholipid"

    # Define a few SMARTS patterns for the phosphocholine headgroup.
    # This pattern looks for a carbon (from glycerol) connected to a phosphate group,
    # which in turn is bound to an oxygen that is linked to an ethylene chain terminating in a trimethylammonium.
    # We include alternative SMARTS to allow for variations in protonation state.
    phosphocholine_smarts_list = [
        "COP(=O)(O)OCC[N+](C)(C)C",      # fully protonated oxygen on phosphate
        "COP(=O)([O-])OCC[N+](C)(C)C",    # one deprotonated oxygen
        "COP([O-])(=O)OCC[N+](C)(C)C"     # alternative order for deprotonation
    ]
    
    # Try to find at least one matching pattern in the molecule.
    headgroup_found = False
    for smarts in phosphocholine_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            headgroup_found = True
            break

    if not headgroup_found:
        return False, "Phosphocholine headgroup substructure not found"

    # If we reach here, we have a phosphorus atom and the phosphocholine headgroup was found.
    # Note: Additional checks (e.g., verifying the glycerol backbone or acyl chain patterns) could
    # be implemented as needed.
    return True, "Molecule contains a phosphocholine moiety (glycerol phosphate ester of phosphocholine)"
    
# Example usage:
if __name__ == "__main__":
    # Test examples (one or two of the provided SMILES strings)
    test_smiles = "P(OC[C@@H](CO)OC(CCCCCCC/C=C\\CC(CCCCCC)O)=O)(=O)(OCC[N+](C)(C)C)[O-]"
    result, reason = is_glycerophosphocholine(test_smiles)
    print(result, reason)