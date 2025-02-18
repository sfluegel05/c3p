"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    Criteria:
      - Must contain a terminal carboxylic acid group.
      - Must have a sufficiently long aliphatic chain (at least 12 carbons).
      - Must contain exactly one epoxide ring (a 3-membered ring with 1 oxygen and 2 carbons).
      - Should not contain additional cyclic structures.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as an epoxy fatty acid, False otherwise.
        str: A reason message explaining the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group.
    # This SMARTS matches a carbon bound to a double-bonded oxygen and an -OH.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found to signify a fatty acid"
    
    # Verify the carboxylic acid is terminal.
    # For a proper fatty acid, the acid carbon should be connected to only one other carbon.
    terminal_found = False
    for match in acid_matches:
        # match[0] is the acid carbon (from C(=O)[O;H])
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Get neighboring atoms that are carbon (ignore oxygen neighbors of the COOH)
        carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_found = True
            break
    if not terminal_found:
        return False, "Carboxylic acid group is not terminal"
    
    # Check that the molecule has a sufficient number of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, f"Too few carbon atoms ({c_count}) to be a typical fatty acid"
    
    # (Optional) Check molecular weight. Most fatty acids weigh over 200 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.2f} Da) too low for a fatty acid"
    
    # Identify rings in the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    epoxide_count = 0
    # Iterate over each ring and look for 3-membered rings (size 3) with exactly one oxygen and two carbons.
    for ring in ring_info:
        if len(ring) == 3:
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            num_oxygen = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
            num_carbon = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            if num_oxygen == 1 and num_carbon == 2:
                epoxide_count += 1
    
    if epoxide_count != 1:
        return False, f"Expected exactly one epoxide ring, found {epoxide_count}"
    
    # To avoid false positives, ensure there are no additional rings.
    # In an epoxy fatty acid, aside from the one epoxide ring, the structure should be largely acyclic.
    total_rings = len(ring_info)
    if total_rings > 1:
        return False, f"Found additional rings aside from the epoxide ring (total rings: {total_rings})"
    
    return True, "Contains a terminal carboxylic acid group, a long aliphatic chain, and one epoxide ring indicative of an epoxy fatty acid"

# Example usage (uncomment to test):
# test_smiles = "C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O"  # (5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid
# result, reason = is_epoxy_fatty_acid(test_smiles)
# print(result, reason)